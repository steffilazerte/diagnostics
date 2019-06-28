#' QQ normality plots with ggplot2
#'
#' Create QQ norm plots through ggplot2
#'
#' @import ggplot2
#' @export
ggQQ <- function(model, level = "all", plot = TRUE, title = NULL) {

  if(!is_lm(model) & !is_lme(model) & !is_mer(model) & !is_gls(model)) {
    stop("'model' must be either a lm, lme, lmer, or glmer model")
  }

  if((is_lm(model) | is_gls(model)) & level != "all") {
    message("levels only apply to mixed models, reverting to level = \"all\"")
    level <- "all"
  }

  # Get residuals in a data frame
  if(level == "all") {
    y <- residuals(model)[!is.na(residuals(model))]
    m <- data.frame(n=1, resid = residuals(model))
    if(is.null(title)) title <- ifelse(is_lm(model) | is_gls(model), "QQ Plot", "QQ Plot: Fixed")
  } else {
    if(is_lme(model)) rand <- nlme::ranef(model, level = level)[[1]]
    if(is_mer(model)) {
      rand <- lme4::ranef(model)
      rand <- rand[[names(rand)[grepl(paste0("^", level, "$"), names(rand))]]][[1]]
    }
    y <- rand[!is.na(rand)]
    m <- data.frame(n = 1, resid = rand)
    if(is.null(title)) title <- paste0("QQ Plot: Random - ", level)
  }

  if(nrow(m) < 2) {
    message("Fewer than 2 levels, no plot possible")
  } else if(sd(m$resid) == 0) {
    message(paste0("Too little variability in '", level, "': QQ Norm plot skipped"))
  } else {
    # Setup for qqline
    x <- qnorm(c(0.25, 0.75))
    y <- quantile(y, c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]

    # Shapiro test statistic
    if(nrow(m) <= 2) {
      message(paste0("Less than 3 data points in '", level, "': Can't compute Shapiro Test"))
      s <- data.frame(p = NA, W = NA)
    } else {
      s <- data.frame(p = round(shapiro.test(m$resid)$p.value, digits = 3),
                      W = round(shapiro.test(m$resid)$statistic, digits = 3))
    }


    if(plot == TRUE){
      p <- ggplot(data = m, aes(sample = resid)) +
        theme_bw() +
        stat_qq(alpha = 0.5) +
        geom_abline(slope = slope, intercept = int, color="blue") +
        labs(y = "Sample Residuals",
             x = "Theoretical Quantiles",
             title = title) +
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



#' Residual plots with ggplot2
#'
#' Create residual plots through ggplot2
#'
#' @import ggplot2
#' @export
ggResid <- function(model) {

    model <- data.frame(fitted = fitted(model),
                        residuals = scale(residuals(model), center=F),
                        row.names = NULL)

    g <- ggplot(model, aes(x = fitted, y = residuals)) +
      theme_bw() +
      geom_point() +
      geom_hline(yintercept=0) +
      labs(title = "Residual Plot")

    return(g)
}

#' @export
diag_plots <- function(model, title = NULL) {

  if(!is_lm(model) & !is_lme(model) & !is_mer(model) & !is_gls(model)) stop("'model' must be either a lm, lme, lmer, or glmer model")

  # Get residual plots and normality plots for all levels of the random variables (if present)
  g <- list()

  # Fixed Effects
  if(!is.null(title)) t_resid <- paste(title, "- Residual Plot")
  g[[1]] <- ggResid(model) + labs(title = t_resid) # Residual
  g[[2]] <- ggQQ(model) + labs(title = "QQ Normality Plot") # Normality

  # Random Effects Normality
  if(is_lme(model) | is_mer(model)) {
    for(a in get_random(model)) g[[length(g) + 1]] <- ggQQ(model, level = a)
  }

  # Show all together
  return(do.call(gridExtra::grid.arrange, c(g, nrow = 1)))  ## Plot
}
