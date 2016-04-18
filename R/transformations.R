
#' Transformations
#'
#' Function to easily transform variables in a variety of ways.
#' By default, variables are anchored at one to maximize the effectiveness of the transformations.
#'
#' @export
trans <- function(x, trans, boxcox = NA, anchor = TRUE){
  if(trans == "boxcox" & is.na(boxcox)) stop("You need to specify a boxcox transformation (lambda)")
  if(length(trans) > 1) stop("You can only specify one transformation at a time")
  if(trans == "boxcox" & !is.na(boxcox)){
    if(boxcox == 0) {
      message("Boxcox of 0, returning log transformation")
      boxcox <- NA
      trans <- "log"
    }
    if(boxcox == 1) {
      message("Boxcox of 1, means no transformation")
      trans = "none"
    }
  }

  if(anchor == TRUE & !(trans %in% c("asn", "none")))  x <- x - min(x, na.rm = TRUE) + 1
  if(trans == "boxcox") return((x ^ boxcox - 1) / boxcox)
  if(trans == "log") return(log(x))
  if(trans == "asn") return(asin(sqrt(x)))
  if(grepl("\\^", trans)) return((x)^eval(parse(text = gsub("\\^", "", trans))))
  if(trans == "inverse") return(1 / x)
  if(trans == "none") return(x)
}

#' Back transform
#'
#' @export
trans.back <- function(x, trans, min = NA, boxcox = NA) {
  if(trans == "log") x <- exp(x)
  if(trans == "asn") x <- sin(x) ^ 2
  if(grepl("\\^", trans)) x <- (x) ^ (1 / eval(parse(text = gsub("\\^", "", trans))))
  if(trans == "boxcox" & !is.na(boxcox)) x <- (x * boxcox + 1) ^ (1 / boxcox)
  if(trans == "inverse") x <- 1 / x
  if(trans == "none") x <- x

  if(!is.na(min) | trans == "asn") x <- x + min - 1
  return(x)
}

#' Plot possible transformations
#'
#' @import gridExtra
#' @export
trans.plot <- function(model, t = c("log","inverse","^1/2","^2","^3"), boxcox = NA, find.best = FALSE, verbose = FALSE) {

  y <- as.character(terms(model)[[2]])

  if(class(model) == "lm") {
    data <- model$model
  } else if (grepl("mer", class(model))) {
    data <- model@frame
  } else {
    data <- nlme::getData(model)
  }

  y.data <- list()
  y.data[[1]] <- data[, y]

  if(find.best == TRUE) {
    if(verbose) message("Finding best Box-Cox...")
    # First Pass
    b <- c(seq(-10, 10, 2))
    b <- b[b != 0]
    best <- data.frame()

    if(verbose) message("\tStarting round 1...")
    for(i in b) {
      temp <- data
      temp[, y] <- trans(y.data[[1]], trans = "boxcox", boxcox = i)
      W <- ggQQ(my.update(model, data = temp), plot = FALSE)$W
      best <- rbind(best, data.frame(b = i, W = W))
      if(verbose) message(paste0("\t\tBoxcox: ", i, "; W = ", W))
    }
    best <- best[order(best$W, decreasing = TRUE), ]
    if(verbose) message(paste0("\t\tBest is: ", best$b[1]))

    # Narrow it down, at most 10 passes
    a <- 1
    while(((best$W[1] - best$W[2]) > 0.001) & a < 10) {
      if(verbose) message(paste0("\tStarting round ", a + 1, "..."))
      n <- c(max(best$b[best$b < best$b[1]]), min(best$b[best$b > best$b[1]][1]))
      if(any(is.na(n))) n <- c(best$b[1] - 0.05, best$b[1] + 0.05)

      best <- data.frame()
      b <- seq(min(n), max(n), length.out = 5)
      b <- b[b != 1]
      for(i in b)  {
        temp <- data
        temp[, y] <- trans(y.data[[1]], trans = "boxcox", boxcox = i)
        W <- ggQQ(my.update(model, data = temp), plot = FALSE)$W
        best <- rbind(best, data.frame(b = i, W = W))
        if(verbose) message(paste0("\t\tBoxcox: ", i, "; W = ", W))
      }
      best <- best[order(best$W, decreasing = TRUE), ]
      if(verbose) message(paste0("\t\tBest is: ", best$b[1]))
      a <- a + 1
    }

    message(paste0("Best Box-cox is: ", best$b[1]))

    t <- "boxcox"
    boxcox <- sort(best$b[1:5])
  }

  if(any(t == "boxcox") & any(is.na(boxcox))) boxcox <- c(-3, -2, -1, 0.5, 2, 3)
  if(any(!is.na(boxcox))) t <- c(t, paste0("Box-cox: ", boxcox))
  t <- t[t != "boxcox"]

  if(verbose) message(paste0("Trying transformations: ", paste0(t, collapse = ", ")))
  for(i in 1:length(t)){
    if(t[i] == "asn") {
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = t[i])
    } else if(!grepl("Box-cox: ", t[i])) {
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = t[i])
    } else {
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = "boxcox", boxcox = as.numeric(gsub("Box-cox: ", "", t[i])))
    }
  }

  q <- list()
  title <- c("No trans", t)
  for(i in 1:length(y.data)) {
    temp <- data
    temp[, y] <- y.data[[i]]
    q[[i]] <- ggQQ(my.update(model, data = temp), title = title[i])
  }

  do.call("grid.arrange", args = c(q, nrow = 1))
}
