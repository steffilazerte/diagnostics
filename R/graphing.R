#' QQ normality plots with ggplot2
#'
#' @import ggplot2
#' @export
ggQQ <- function(model, level = "all", plot = TRUE, title = NULL) {

  cls <- class(model)
  if(grepl("mer", cls)) cls <- "lmer"

  if(!(cls %in% c("lm", "lme", "lmer"))) {
    stop("'model' must be either a lm, lme, lmer, or glmer model")
  }

  if(cls == "lm" & level != "all") {
    message("levels only apply to mixed models, reverting to level = \"all\"")
    level <- "all"
  }

  # Get residuals in a data frame
  if(level == "all") {
    y <- resid(model)[!is.na(resid(model))]
    m <- data.frame(n=1, resid = resid(model))
    title <- ifelse(cls == "lm", "QQ Plot", "QQ Plot: Fixed")
  } else {
    if(cls == "lme") rand <- nlme::ranef(model, level = level)[[1]]
    if(cls == "lmer") {
      rand <- lme4::ranef(model)
      rand <- rand[[names(rand)[grepl(paste0("^", level, "$"), names(rand))]]][[1]]
    }
    y <- rand[!is.na(rand)]
    m <- data.frame(n = 1, resid = rand)
    title <- paste0("QQ Plot: Random - ", level)
  }

  if(nrow(m) < 2) {
    message("Fewer than 2 levels, no plot possible")
  } else {
    # Setup for qqline
    x <- qnorm(c(0.25,0.75))
    y <- quantile(y, c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]

    # Shapiro test statistic
    if(nrow(m) > 2 | all(m$resid == 0)) {
      s <- data.frame(p = round(shapiro.test(m$resid)$p.value, digits = 3),
                      W = round(shapiro.test(m$resid)$statistic, digits = 3))
    } else if(nrow(m) <= 2) {
      message(paste0("Less than 3 data points in '", level, "': Can't compute Shapiro Test"))
      s <- data.frame(p = NA, W = NA)
    } else if(sd(m$resid) == 0) {
      message(paste0("Too little variability in '", level, "': Can't compute Shapiro Test"))
      s <- data.frame(p = NA, W = NA)
    }


    if(plot == TRUE){
      p <- ggplot(data = m, aes(sample = resid)) +
        theme_bw() +
        stat_qq(alpha = 0.5) +
        geom_abline(slope = slope, intercept = int, color="blue") +
        labs(y = "Sample Residuals",
             x = "Theoretical Quantiles",
             title = ifelse(is.null(title), "QQ Plot", title)) +
        annotate(geom = "text",
                 x = -Inf, y = +Inf,
                 hjust = 0, vjust = 1,
                 size = 3,
                 label = paste0("Shapiro Test \np = ", s$p, "\nW = ", s$W))
      return(p)
    } else {
      return(s)
    }
  }
}



#' QQ residual plots with ggplot2
#'
#' @import ggplot2
#' @export
ggResid <- function(model) {

    model <- data.frame(fitted = fitted(model),
                        residuals = scale(resid(model), center=F),
                        row.names = NULL)

    g <- ggplot(model, aes(x = fitted, y = residuals)) +
      theme_bw() +
      geom_point() +
      geom_hline(yintercept=0) +
      labs(title = "Residual Plot")

    return(g)
}
