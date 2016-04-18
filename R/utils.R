# Hadley's Model Update function
my.update <- function(mod, formula = NULL, data = NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  eval(call, env, parent.frame())
}

test.plot <- function(model, level = "all", type = "QQ") {
  par(mfrow=c(1, 2), mar = c(2,2,2,2))
  if(type == "QQ") {
    if(level == "all"){
      qqnorm(resid(model), main = "QQ Plot")
      qqline(resid(model))
    } else {
      qqnorm(ranef(model)[[level]][[1]], main = "QQ Plot")
      qqline(ranef(model)[[level]][[1]])
    }
  } else if (type == "R") plot(fitted(model), residuals(model))
  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
  vp1 <-plotViewport(c(0,0,0,0))
  if(type == "QQ") g <- ggQQ(model, level = level) else if(type == "R") g <- ggResid(model)
  print(g, vp = vp1)
}
