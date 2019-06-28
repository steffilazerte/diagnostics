# From Hadley's Model Update function
# http://stackoverflow.com/a/13690928/3362144
model_update <- function(model, formula = NULL, data = NULL) {
  call <- stats::getCall(model)
  if (is.null(call)) stop("Model object has no call", call. = FALSE)
  term <- stats::terms(model)
  if (is.null(term)) stop("Model object has no terms", call. = FALSE)

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- stats::update.formula(call$formula, formula)

  env <- attr(term, ".Environment")
  e <- eval(call, env, parent.frame())

  # Convert back to lmerTest if was to begin
  if(any(grepl("Test", class(model)))) e <- lmerTest::as_lmerModLmerTest(e)

  return(e)
}

# Fixed effects
get_fixed <- function(model){
 if(is_lm(model)) return(stats::coef(model))
 if(is_lme(model)) return(nlme::fixef(model))
 if(is_mer(model)) return(lme4::fixef(model))
}

# Retrieve data from model objects
get_data <- function(model){
  call <- stats::getCall(model)
  term <- stats::terms(model)
  env <- attr(term, ".Environment")
  return(eval(call$data, envir = env))
}

is_lm <- function(model) if(any(class(model) %in% c("lm", "glm"))) return(TRUE) else return(FALSE)
is_gls <- function(model) if(any(class(model) %in% "gls")) return(TRUE) else return(FALSE)
is_lme <- function(model) if(any(class(model) %in% "lme")) return(TRUE) else return(FALSE)
is_mer <- function(model) if(any(grepl("merMod", class(model)))) return(TRUE) else return(FALSE)
is_merTest <- function(model) if(any(grepl("LmerTest", class(model)))) return(TRUE) else return(FALSE)

get_random <- function(model) {
  if(is_lme(model)) return(names(model$modelStruct[[1]]))
  if(is_mer(model)) return(names(ranef(model)))
}

get_ncoef <- function(model) {
  if(is_lm(model)) return(length(model$coefficients))
  if(is_lme(model)) return(length(model$coefficients$fixed))
  if(is_mer(model)) return(ncol(summary(model)$vcov))
}

# Get difference in p-values bins
diff_p <- function(new, orig, model) {
  new %>%
    model_update(model = model, data = .) %>%
    get_table() %>%
    dplyr::mutate(diff5 = (P < 0.05) != (orig$P < 0.05),
                  diff10 = (P < 0.10) != (orig$P < 0.10))
}
