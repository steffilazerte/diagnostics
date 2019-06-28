#' Universal calculation of Variance Inflation Factors (VIF)
#'
#' Calculates VIF for most models
#'
#' For non-mixed models (lm and glm) uses vif from car package. For mixed models
#' uses function from Austin F. Frank:
#' https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
#'
#' "VIF values greater than ten suggest strong collinearity." (Quinn & Keough, 2002)
#'
#' @export
uni_vif <- function(model) {
  if(is_lm(model)) return(car::vif(model))
  if(is_lme(model) | is_mer(model)) return(vif_mer(model))
}

# Based on function by Austin F. Frank:
# https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
vif_mer <- function(model) {
    ## adapted from rms::vif
    v <- as.matrix(stats::vcov(model))
    nam <- names(get_fixed(model))

    # exclude intercepts
    ns <- which(nam == "(Intercept)")
    if (ns > 0) {
       v <- v[-(1:ns), -(1:ns), drop = FALSE]
       nam <- nam[-(1:ns)]
    }

    d <- diag(v)^0.5
    v <- Matrix::diag(solve(v/(d %o% d)))
    names(v) <- nam
    return(v)
}

#' Universal calculation of Condition Number (kappa)
#'
#' Calculates Condition Number of matrix for common model types
#'
#' For non-mixed models (lm and glm) uses kappa() from base package. For mixed models
#' uses function from Austin F. Frank:
#' https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
#'
#' Condition numbers should be less than 10 or 20 or you have multicollinearity
#' "values greater than 30 indicate collinearities that require attention." (Quinn and Keough
#' 2002)
#'
#' Note that kappa(lm) == uni_kappa(lm) but kappa(glm) != uni_kappa(glm)
#' Has to do with method of calculation, use with caution
#'
#' @export
uni_kappa <- function(model, scale = FALSE, center = FALSE, exact = FALSE, add_intercept = TRUE) {

  # Get Model Matrices
  if(is_lm(model) | is_lme(model)) X <- stats::model.matrix(model)
  if(is_mer(model)) X <- lme4::getME(model, name = "X")

  # Keep names
  nam <- names(get_fixed(model))

  # Exclude intercepts
  intercept <- which(nam == "(Intercept)")
  if(intercept > 0) {
    X <- X[, -(1:intercept), drop = FALSE]
    nam <- nam[-(1:intercept)]
  }

  if(add_intercept) {
    X <- cbind(rep(1), scale(X, scale = scale, center = center))
    kappa(X, exact = exact)
  } else {
    kappa(scale(X, scale = scale, center = scale), exact = exact)
  }
}

#' Overdispersion for glmer
#'
#' Approximate test of overdispersion from:
#' http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
#'
#' @author Ben Bolker
#'
#' @export
overdisp <- function(model) {
  ## number of variance parameters in
  ## an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(lme4::VarCorr(model), vpars)) + length(get_fixed(model))
  rdf <- nrow(stats::model.frame(model)) - model.df
  rp <- stats::residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- stats::pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  return(c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval))
}

#' Test for influential observations
#'
#'
#' @export
sig_test <- function(model, group = "obs", all = FALSE, verbose = FALSE){

  if(!is_lm(model) & !is_lme(model) & !is_mer(model)) stop("'model' must be either a lm, glm, lme, lmer, or glmer model")

  if(is_mer(model) & !is_merTest(model)){
    stop("Require lmer/glmer models to be run with the lmerTest package to calculate degrees of freedom. Load the lmerTest package with 'library(lmerTest)' and run the lmer()/glmer() model again.")
  }

  data <- as.data.frame(get_data(model))
  orig <- get_table(model)

  # Check that obs correct
  if(is_lm(model) & group != "obs") {
    message("Different groupings only apply to mixed models, reverting to group = \"obs\"")
    group <- "obs"
  } else if (group != "obs" & !(group %in% names(data))) {
    stop("'group' must be the name of a column (i.e. 'ID'), or 'obs', to reflect observation level grouping.")
  }

  # Group and omit group observations from lists of data sets
  if(group == "obs") {
    g <- tibble::tibble(type = "obs", obs = 1:nrow(data)) %>%
      dplyr::group_by(obs) %>%
      dplyr::do(data = data[-.$obs, ]) %>%
      dplyr::ungroup()
  } else {
    g <- tibble::tibble(type = group, obs = unique(data[, group])) %>%
      dplyr::group_by(obs) %>%
      dplyr::do(data = data[data[, group] != .$obs, ]) %>%
      dplyr::ungroup()
  }

  # Get new sig with sequential omission
  tests <- g %>%
    dplyr::mutate(new_fit = purrr::map(data, diff_p, orig = orig, model = model)) %>%
    tidyr::unnest(new_fit)

  if(any(tests$P == "")) {
    message("Influential: Some models had problems and alternate P-values could not be computed. Influential tests unreliable...")
  } else {

    if(!all) tests <- dplyr::filter(tests, diff5 | diff10)

    tests <- tests %>%
      dplyr::left_join(orig, by = "Parameter", suffix = c(".new", ".orig")) %>%
      dplyr::select(Parameter, Value.new, SE.new, P.new, Value.orig, SE.orig, P.orig)

    if(nrow(tests) == 0) message("No influential observations") else {
      return(tests)
      if(group == "obs") message("Note that for group = \"obs\" numbers returned refer to the row in a data frame without NA values.")
    }
  }
}

#' Test for multicolinearity
#'
#' Test for MC
#'
#' @export
multicol <- function(model) {

  if(get_ncoef(model) > 2){
    if(is_lm(model)) {
      v <- uni_vif(model)
      k <- uni_kappa(model)
    } else {
      v <- uni_vif(model)
      k <- uni_kappa(model)
    }
  } else {
    v <- k <- NA
    message("Can't calculate multicolinearity if there are less than 2 predictors.")
  }
  return(list(v = v, k = k))
}


#' Run Diagnostics
#'
#' Run complete diagnostics
#'
#' @import gridExtra
#' @import influence.ME
#' @import car
#' @export
diagnostic <- function(model, group = "obs", graphs = TRUE, influence = TRUE, multicol = TRUE, verbose = FALSE) {

  if(!is_lm(model) & !is_lme(model) & !is_mer(model)) stop("'model' must be either a lm, lme, lmer, or glmer model")

  # Get residual plots and normality plots for all levels of the random variables (if present)
  if(graphs == T){
    g <- list()

    # Fixed Effects
    g[[1]] <- ggResid(model) + labs(title = "Residual Plot") # Residual
    g[[2]] <- ggQQ(model) + labs(title = "QQ Normality Plot") # Normality

    # Random Effects Normality
    if(is_lme(model) | is_mer(model)) {
      for(a in get_random(model)) g[[length(g) + 1]] <- ggQQ(model, level = a)
    }

    # Show all together
    do.call(gridExtra::grid.arrange, c(g, nrow = 1))  ## Plot
  }

  # Get influential observations
  if(influence == T){
    # Influential obs
    i <- sig_test(model, group = group, verbose = verbose)
    if(!is.null(i)) {
      i <- i %>%
        dplyr::select(Parameter, Value.orig, Value.new, P.orig, P.new) %>%
        dplyr::mutate(Obs = group,
                      Diff = NA,
                      Diff = replace(Diff, P.orig < 0.05 & (P.new >= 0.05 & P.new < 0.10), "Sig => Trend"),
                      Diff = replace(Diff, P.orig < 0.05 & (P.new >= 0.10), "Sig => Non Sig"),
                      Diff = replace(Diff, (P.orig >= 0.05 & P.orig < 0.10) & P.new < 0.05, "Trend => Sig"),
                      Diff = replace(Diff, (P.orig >= 0.05 & P.orig < 0.10) & P.new >= 0.10, "Trend => Non Sig"),
                      Diff = replace(Diff, P.orig >= 0.10 & P.new < 0.05, "Non Sig => Sig"),
                      Diff = replace(Diff, P.orig >= 0.10 & (P.new >= 0.05 & P.new < 0.10), "Non Sig => Trend")) %>%
        dplyr::rename(Param = Parameter, Value = Value.orig, `New Value` = Value.new, P = P.orig, `New P` = `P.new`)

    }

    # Cooks distances
    #cooks <- predictmeans::CookD(model, group = if(group == "obs") NULL else group, plot = FALSE)
    #cooks <- data.frame(cook = cooks, obs = names(cooks))
    #cooks <- cooks[cooks > 1]
  }

  # Get Multicolinearity
  if(multicol == T){
    m <- multicol(model)
    k <- m$k
    v <- m$v
  }

  ## Display
  if(influence == T){
    if(!is.null(i)) {
    message("Influence: \n")
    print(i)
   # message(paste0("Cooks Distances: ", paste0(cooks, collapse = "; ")))
    }
  }
  if(multicol == T) {
    if(any(!is.na(v), !is.na(k))) {
    message("Multicolinearity:\n")
    print(paste0("VIF: ", paste(names(v), round(v, digits = 1), sep = " ")))
    print(paste0("Kappa: ", round(k,digits = 1)))
    }
  }
}

#' Report summary statistics for different models
#'
#' Grab data from differnet model types
#'
#' @export
get_table <- function(model, analysis = NULL, type = "summary", level = NULL, pca = 1){
  if(type == "summary") {
    if(is_lm(model)){
      x <- as.data.frame(summary(model)$coefficients)
      names(x) <- c("Value","SE","T","P")
      x$df <- max(summary(model)$df)
      i1 <- confint(model)
      i2 <- confint(model, level = 0.9)
      x$n <- length(model$fitted.values)
      R2 <- c(summary(model)$r.squared, "")
    } else if(is_lme(model)) {
      x <- as.data.frame(summary(model)$tTable)
      names(x) <- c("Value","SE","df", "T","P")
      i1 <- nlme::intervals(model, which = "fixed")[[1]]
      i2 <- nlme::intervals(model, which = "fixed", level = 0.9)[[1]]
      if(is.null(level)) x$n <- nrow(summary(model)$groups) else x$n <- length(unique(model$data[,level]))
      R2 <- MuMIn::r.squaredGLMM(model)
    } else if(is_mer(model)) {
      if(class(model) == "glmerMod") {
        x <- as.data.frame(summary(model)$coefficients)
        names(x) <- c("Value", "SE", "T","P")
        x$note <- "Test statistic is z not T"
      } else if (grepl("Test", class(model))) {
        x <- as.data.frame(lmerTest:::summary.lmerModLmerTest(model)$coefficients)
        if(any(names(x) == "df")) {
          names(x) <- c("Value", "SE", "df","T","P")
        } else names(x) <- c("Value", "SE", "T")
      } else {
        x <- as.data.frame(summary(model)$coefficients)
        names(x) <- c("Value", "SE", "T")
      }
      R2 <- MuMIn::r.squaredGLMM(model)
      i1 <- data.frame(confint(model, method = "Wald"))
      i1 <- i1[-grep(".sig",row.names(i1)),]
      i2 <- data.frame(confint(model, method = "Wald", level = 0.9))
      i2 <- i2[-grep(".sig",row.names(i2)),]
      if(is.null(level)) x$n <- nrow(model@frame) else x$n <- length(unique(model@frame[,level]))
    }
    x$Parameter <- row.names(x)
    x$CI.95 <- x$Value - i1[,1]
    x$CI.90 <- x$Value - i2[,1]
    x$sig.95 <- ifelse((abs(x$Value) - x$CI.95) > 0, T, F)
    x$sig.90 <- ifelse((abs(x$Value) - x$CI.90) > 0, T, F)
    if(!any(grepl("^df$", names(x)))) x$df <- ""
    if(!any(grepl("^P$", names(x)))) x$P <- ""
    if(!any(grepl("^note$", names(x)))) x$note <- ""
    x$R2.m <- R2[1]
    x$R2.c <- R2[2]
    x <- x[,c("Parameter","Value","CI.95","df","T","P","n","SE","sig.95","CI.90","sig.90","R2.m", "R2.c","note")]
  } else if(type == "anova") {
    if(is_lme(model)) {
      x <- anova(model)
      names(x) <- c("df_num","df_den","F","P")
      x <- cbind(Parameter = row.names(x), x)
    } else if(is_mer(model)) {
      x <- anova(model)
      names(x) <- c("df","SS","MS","F")
      x <- cbind(Parameter = row.names(x), x)
    }
  }
  if(!is.null(analysis)) x <- cbind(Analysis = analysis, x)
  return(x)
}
