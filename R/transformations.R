

trans_basic <- function(x, trans, boxcox = NA) {
  if(trans == "boxcox") x <- (x ^ boxcox - 1) / boxcox
  if(trans == "log") x <- log(x, base = 10)
  if(trans == "ln") x <- log(x)
  if(trans == "asn") x <- asin(sqrt(x))
  if(grepl("\\^", trans)) x <- (x)^eval(parse(text = gsub("\\^", "", trans)))
  if(trans == "inverse") x <- 1 / x
  if(trans == "none") x <- x
  return(x)
}


#' Transformations
#'
#' Function to easily transform variables in a variety of ways.
#' By default, variables are anchored at one to maximize the effectiveness of the transformations.
#'
#' @export
trans <- function(x, trans, boxcox = NA, anchor = TRUE, centre = FALSE){
  if(trans == "boxcox" & is.na(boxcox)) stop("You need to specify a boxcox transformation (lambda)")
  if(length(trans) > 1) stop("You can only specify one transformation at a time")
  if(trans == "asn" & !all(x >= 0 & x <= 1) & centre == FALSE) stop("For Arc-Sin Square-root tranformations, data needs to be proportions(i.e. range from 0 to 1); If you have a change in proportion (i.e. range from -1 to 1) specify centre = TRUE.")
  if(trans == "asn" & !all(x >= (-1) & x <= 1) & centre == TRUE) stop("For Arc-Sin Square-root tranformations, data needs to be proportions: range from 0 to 1; or, if you have a change in proportions: range from -1 to 1")
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

  # If centred, grab which negative then abolute the whole thing
  if(centre == TRUE & !(trans %in% c("inverse", "none"))) {
    neg <- which(x < 0) # Save negative locs
    x <- abs(x) # trans to absolute values
  }

  # Anchor data
  if(anchor == TRUE & !(trans %in% c("asn", "none")))  x <- x - min(x, na.rm = TRUE) + 1  ## anchor at 1

  # Transform
  x <- trans_basic(x, trans = trans, boxcox = boxcox)

  # If centred, add in negatives around zero
  if(centre == TRUE & !(trans %in% c("inverse", "none"))) x[neg] <- -1 * x[neg]
  return(x)
}



#' Back transform
#'
#' @export
trans.back <- function(x, trans, boxcox = NA, min = NA, centre = FALSE) {

  if(centre == TRUE & !(trans %in% c("inverse", "none"))) {
    neg <- which(x < 0)
    x <- abs(x)
    min <- abs(min)
  }

  if(trans == "log") x <- 10^(x)
  if(trans == "ln") x <- exp(x)
  if(trans == "asn" & centre == FALSE) x <- sin(x) ^ 2
  if(trans == "asn" & centre == TRUE) x <- c(-sin(-x[x < 0])^2, sin(x[x >= 0])^2)
  if(grepl("\\^", trans)) x <- (x) ^ (1 / eval(parse(text = gsub("\\^", "", trans))))
  if(trans == "boxcox" & !is.na(boxcox)) x <- (x * boxcox + 1) ^ (1 / boxcox)
  if(trans == "inverse") x <- 1 / x
  if(trans == "none") x <- x

  if(!is.na(min) & !(trans %in% c("asn", "none"))) x <- x + min - 1

  if(centre == TRUE & !(trans %in% c("inverse", "none"))) x[neg] <- -1 * x[neg]

  return(x)
}

#' Plot possible transformations
#'
#' @import gridExtra
#' @export
trans.plot <- function(model, t = c("log","inverse","^1/2","^2","^3"), boxcox = NA, find.best = FALSE, centre = FALSE, verbose = FALSE) {

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
      temp[, y] <- trans(y.data[[1]], trans = "boxcox", boxcox = i, centre = centre)
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
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = t[i], centre = centre)
    } else if(!grepl("Box-cox: ", t[i])) {
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = t[i], centre = centre)
    } else {
      y.data[[(i + 1)]] <- trans(y.data[[1]], trans = "boxcox", boxcox = as.numeric(gsub("Box-cox: ", "", t[i])), centre = centre)
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
