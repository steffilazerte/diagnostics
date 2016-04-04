
#' Transformations
#'
#' @export
trans <- function(x, trans, boxcox=NA, anchor = TRUE){
  if(anchor == T & trans != "asn" & (is.na(boxcox) | boxcox!=1))  x <- x - min(x, na.rm=T) + 1
  if(trans == "log") return(log(x))
  if(trans == "asn") return(asin(sqrt(x)))
  if(grepl("\\^", trans)) return((x)^eval(parse(text=gsub("\\^","",trans))))
  if(trans == "boxcox" & !is.na(boxcox) & boxcox!=1) return((x^boxcox - 1)/boxcox) else return(x)
  if(trans == "inverse") return(1/x)
  if(trans == "none") return(x)
}

#' Back transform
#'
#' @export
trans.back <- function(x, trans, min = NA, boxcox = NA) {
  if(trans == "log") x <- exp(x)
  if(trans == "asn") x <- sin(x)^2
  if(grepl("\\^", trans)) x <- (x)^(1/eval(parse(text=gsub("\\^","",trans))))
  if(trans == "boxcox" & !is.na(boxcox)) x <- (x * boxcox + 1)^(1/boxcox)
  if(trans == "inverse") x <- 1/x
  if(trans == "none") x <- x

  if(!is.na(min) | trans == "asn") x <- x + min - 1
  return(x)
}

#' Plot possible transformations
#'
#' @import ggplot2
#' @import gridExtra
#' @import nlme
#' @import lme4
#' @export
trans.plot <- function(model, t = c("log","inverse","^1/2","^2","^3"), boxcox = NA, find.best = F) {

  y <- as.character(terms(model)[[2]])

  if(class(model) == "lm") {
    data <- model$model
  } else if (grepl("mer", class(model))) {
    data <- model@frame
  } else {
    data <- getData(model)
  }

  y.data <- list()
  y.data[[1]] <- data[,y]

  if(find.best == T) {
    # First Pass
    b <- c(seq(-10, 10, 2))
    b <- b[b!=0]
    b <- sort(b)
    best <- data.frame()
    for(i in b){
      temp <- data
      temp[,y] <- trans(y.data[[1]], trans = "boxcox", boxcox = i)
      best <- rbind(best, data.frame(b = i, W = ggQQ(update(model, data = temp), plot = F)$w))
    }
    best <- best[order(best$W, decreasing = T),]

    # Narrow it down, at most 10 passes
    a <- 1
    while(((best$W[1] - best$W[2]) > 0.001) & a < 10) {
      n <- c(best$b[best$b < best$b[1]][1], best$b[best$b > best$b[1]][1])
      if(any(is.na(n))) n <- c(best$b[1] - 0.05, best$b[1] + 0.05)
      best <- data.frame(b = seq(min(n),max(n), length.out = 5), W = NA)
      for(i in best$b)  {
        temp <- data
        temp[,y] <- trans(y.data[[1]],trans = "boxcox", boxcox = i)
        best$W[best$b == i] <- ggQQ(update(model, data = temp), plot = F)$w
      }
      best <- best[order(best$W, decreasing = T),]
      a <- a + 1
    }
    t <- "boxcox"
    boxcox <- sort(best$b[1:5])
  }

  if(any(t == "boxcox") & any(is.na(boxcox))) boxcox <- c(-3,-2,-1,1,2,3)
  if(any(!is.na(boxcox))) t <- c(t,paste0("boxcox.",boxcox))
  t <- t[t!="boxcox"]

  for(i in 1:length(t)){
    if(t[i] == "asn") {
      y.data[[(i+1)]] <- trans(y.data[[1]], trans = t[i])
    } else if(!grepl("boxcox",t[i])) {
      y.data[[(i+1)]] <- trans(y.data[[1]], trans = t[i])
    } else {
      y.data[[(i+1)]] <- trans(y.data[[1]], trans = "boxcox", boxcox = as.numeric(gsub("boxcox.","",t[i])))
    }
  }

  q <- list()
  t <- gsub("boxcox.", "Box-cox: ", t)
  title <- c("No trans",t)
  for(i in 1:length(y.data)) {
    temp <- data
    temp[,y] <- y.data[[i]]
    q[[i]] <- ggQQ(update(model, data = temp), title = title[i])
  }
  do.call("grid.arrange", args = c(q, nrow = 1))
}
