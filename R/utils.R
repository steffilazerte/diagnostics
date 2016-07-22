# Hadley's Model Update function
my.update <- function(model, formula = NULL, data = NULL) {
  call <- getCall(model)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(model)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  e <- eval(call, env, parent.frame())

 # if(grepl("LmerTest", class(model))) {
#    e <- as(e, "merModLmerTest")
#      #call <- as.list(call)
#    #call[[1]] <- as.symbol("lmer")
#    #call <- as.call(call)
#  }

  return(e)
}

my.lmer <- lmerTest::lmer

#' Get data from a merMod obejct
#'
#' New method
#'
#' @importFrom nlme getData
#' @export
getData.merMod <- function(model){
  call <- getCall(model)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(model)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  env <- attr(term, ".Environment")
  data <- get(as.character(call$data), envir = env)
  return(data)
}

#' Get data from a lm obejct
#'
#' New method
#'
#' @export
getData.lm <- function(model){
  data <- model$model
  return(data)
}

get.random <- function(model) {
  if(class(model) == "lme") return(names(model$modelStruct[[1]]))
  if(grepl("merMod", class(model))) return(names(ranef(model)))
}

get.ncoef <- function(model) {
  if(class(model) == "lm") return(length(model$coefficients))
  if(class(model) == "lme") return(length(model$coefficients$fixed))
  if(grepl("merMod", class(model))) return(ncol(summary(model)$vcov))
}
