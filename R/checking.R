#' KMO for PCA
#'
#' Test KMO of pca analyses
#'
#' @import MASS
#' @export
kmo = function( data ){
  ## FROM http://www.opensubscriber.com/message/r-help@stat.math.ethz.ch/7315408.html
  X <- cor(as.matrix(data))
  iX <- ginv(X)
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a)
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy

  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the
                                        # correlation matrix. That is the
                                        # negative of the partial correlations,
                                        # partialling out all other variables.

  kmo <- BB/(AA+BB)                     # overall KMO statistic

  # Reporting the conclusion
    if (kmo >= 0.00 && kmo < 0.50){
      test <- 'The KMO test yields a degree of common variance
unacceptable for FA.'
    } else if (kmo >= 0.50 && kmo < 0.60){
      test <- 'The KMO test yields a degree of common variance miserable.'
    } else if (kmo >= 0.60 && kmo < 0.70){
      test <- 'The KMO test yields a degree of common variance mediocre.'
    } else if (kmo >= 0.70 && kmo < 0.80){
      test <- 'The KMO test yields a degree of common variance middling.'
    } else if (kmo >= 0.80 && kmo < 0.90){
      test <- 'The KMO test yields a degree of common variance meritorious.'
    } else {
      test <- 'The KMO test yields a degree of common variance marvelous.'
    }

    ans <- list(  overall = kmo,
                  report = test,
                  individual = MSA,
                  AIS = AIS,
                  AIR = AIR )
    return(ans)

}    # end of kmo()


#' Bartlett sphere
#'
#' Bartlett sphere test for PCA
#'
#' @export
bartlett.sphere<-function(data){
  ##FROM https://stat.ethz.ch/pipermail/r-help/2011-June/281243.html
  chi.square <- -( (nrow(data)-1) - (2*ncol(data)-5)/6 )*log(det(cor(data,use='pairwise.complete.obs')))
  paste0('chi.square value ',chi.square , ' on ', (ncol(data)^2-ncol(data))/2, ' degrees of freedom.', ' p-value: ', pchisq(chi.square,(ncol(data)^2-ncol(data))/2,lower.tail=F))
}

#' Check VIF
#'
#' https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
#'
#' @import lme4
#' @import Matrix
#' @export
vifmer <- function (fit) {
    ## adapted from rms::vif

    v <- vcov(fit)
    nam <- names(fixef(fit))

    # exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
       v <- v[-(1:ns), -(1:ns), drop = FALSE]
       nam <- nam[-(1:ns)]
    }

    d <- diag(v)^0.5
    v <- Matrix::diag(solve(v/(d %o% d)))
    names(v) <- nam
    return(v)
}

#' Check Kappa
#'
#' https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
#'
#' Condition number test less than 10 or 20 (less than 30)
#'
#' @import lme4
#' @export
kappamer <- function (fit,
                      scale = TRUE, center = FALSE,
                      add.intercept = FALSE,
                      exact = FALSE) {
  cls <- class(fit)
  if(grepl("merMod", cls)) cls <- "lmer"

  # Get Model Matrices
  if(cls == "lmer") X <- getME(fit, name = "X")
  if(cls == "lme") X <- model.matrix(fit, data = getData(fit))

  nam <- names(fixef(fit))

  # exclude intercepts
  nrp <- sum(1 * (nam == "(Intercept)"))
  if (nrp > 0) {
      X <- X[, -(1:nrp), drop = FALSE]
      nam <- nam[-(1:nrp)]
  }

  if (add.intercept) {
      X <- cbind(rep(1), scale(X, scale = scale, center = center))
      kappa(X, exact = exact)
  } else {
      kappa(scale(X, scale = scale, center = scale), exact = exact)
  }
}

#' Overdispersion glmer
#'
#' Tests for overdispersion
#'
#' @export
overdisp.glmer <- function(modelglmer) {
  ## computing  estimated scale  ( binomial model)
  ##following  D. Bates :
  ##That quantity is the square root of the penalized residual sum of
  ##squares divided by n, the number of observations, evaluated as:
  n <- length(resid(modelglmer))
  return(  sqrt( sum(c(resid(modelglmer), modelglmer@u) ^2) / n ) )
}

#' Overdispersion glmer 2
#'
#' Tests for overdispersion
#'
#' @export
overdisp2.glmer <- function(model) {
  # number of variance parameters in
  #   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#' Test for influential observtions
#'
#' Test for influential obserations
#'
#' @export
sig.test <- function(model, group = "obs", all = FALSE, verbose = FALSE){

  cls <- class(model)
  if(grepl("merMod", cls)) cls <- "lmer"
  if(!(cls %in% c("lm", "lme", "lmer"))) stop("'model' must be either a lm, lme, lmer, or glmer model")

  if(cls == "lmer" & !grepl("LmerTest", class(model)) & !grepl("glmer", class(model))){
    stop("Can't test for influential observations without p-values, right now. Try again with library(lmerTest)")
  }

  data <- getData(model)
  orig <- get.table(model)

  # Check that obs correct
  if(cls == "lm" & group != "obs") {
    message("Different groupings only apply to mixed models, reverting to group = \"obs\"")
    group <- "obs"
  } else if (group != "obs" & !(group %in% names(data))) {
    stop("'group' must be the name of a column (i.e. 'ID'), or 'obs', to reflect observation level grouping.")
  }

  # Grouping
  if(group == "obs") {
    g <- as.list(1:nrow(data))
    names(g) <- 1:nrow(data)
  } else {
    g <- lapply(unique(data[, group]), FUN = function(x) which(data[, group] == x))
    names(g) <- unique(data[, group])
  }

  # Get new sig with sequential omission
  if(cls != "lmer") tests <- plyr::ldply(g, .fun = function(x) get.table(my.update(model, data = data[-x, ]))[, c("Parameter", "Value", "P")], .id = 'group', .progress = ifelse(verbose, "text", "none"))
  if(cls == "lmer") tests <- plyr::ldply(g, .fun = function(x) get.table(methods::as(my.update(model, data = data[-x, ]), "merModLmerTest"))[, c("Parameter", "Value", "P")], .id = 'group', .progress = ifelse(verbose, "text", "none"))

  if(any(tests$P == "")) {
    message("Influential: Some models had problems and alternate P-values could not be computed. Influential tests unreliable...")
  } else {
    output <- merge(tests, orig[, c("Parameter", "Value", "P")], by = "Parameter", suffixes = c("", ".orig"))
    output[, c("Value", "P", "Value.orig", "P.orig")] <- apply(output[, c("Value", "P", "Value.orig", "P.orig")], 2, round, 4)
    output$diff5 <- (output$P < 0.05) != (output$P.orig < 0.05)
    output$diff10 <- (output$P < 0.10) != (output$P.orig < 0.10)

    if(!all) output <- output[output$diff5 | output$diff10, ]

    if(nrow(output) == 0) message("No influential observations") else {
      return(output)
      message("Note that for group = \"obs\" numbers returned refer to the row in a data frame without NA values.")
    }
  }
}

#' Test for multicolinearity
#'
#' Test for MC
#'
#' @export
multicol <- function(model) {
  cls <- class(model)
  if(grepl("merMod", cls)) cls <- "lmer"

  if(get.ncoef(model) > 2){
    if(cls == "lm") {
      v <- car::vif(model)[, 1]
      k <- kappa(model)
    }
    if(cls != "lm") {
      v <- vifmer(model)
      k <- kappamer(model)
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

  cls <- class(model)
  if(grepl("merMod", cls)) cls <- "lmer"
  if(!(cls %in% c("lm", "lme", "lmer"))) stop("'model' must be either a lm, lme, lmer, or glmer model")

  #require(nlme)
  #require(nlmeU)

  ## Get residual plots and normality plots for all levels of the random variables (if present)
  if(graphs == T){
    g <- list()

    # Fixed Effects
    g[[1]] <- ggResid(model) + labs(title = "Residual Plot") # Residual
    g[[2]] <- ggQQ(model) + labs(title = "QQ Normality Plot") # Normality

    # Random Effects Normality
    if(cls %in% c("lme", "lmer")) {
      for(a in get.random(model)) g[[length(g) + 1]] <- ggQQ(model, level = a)
    }

    # Show all together
    do.call(gridExtra::grid.arrange, c(g, nrow = 1))  ## Plot
  }

  # Get influential observations
  if(influence == T){
    # Influential obs
    i <- sig.test(model, group = group, verbose = verbose)
    if(!is.null(i)) {
      i <- i[, c("group", "Parameter", "Value.orig", "Value", "P.orig", "P")]
      i$Diff <- NA
      i$Diff[i$P.orig < 0.05 & (i$P >= 0.05 & i$P < 0.10)] <- "Sig => Trend"
      i$Diff[i$P.orig < 0.05 & (i$P >= 0.10)] <- "Sig => Non Sig"
      i$Diff[(i$P.orig >= 0.05 & i$P.orig < 0.10) & (i$P < 0.05)] <- "Trend => Sig"
      i$Diff[(i$P.orig >= 0.05 & i$P.orig < 0.10) & (i$P >= 0.10)] <- "Trend => Non Sig"
      i$Diff[(i$P.orig >= 0.10) & (i$P < 0.05)] <- "Non Sig => Sig"
      i$Diff[(i$P.orig >= 0.10) & (i$P >= 0.05 & i$P < 0.10)] <- "Non Sig => Trend"
      names(i) <- c("Obs", "Param", "Value", "New Value", "P", "New P", "Diff")

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
get.table <- function(model, analysis = NULL, type = "summary", level = NULL, pca = 1){
  if(type == "summary") {
    cls <- class(model)
    if(grepl("merMod", cls)) cls <- "lmer"
    if(cls == "lm"){
      x <- as.data.frame(summary(model)$coefficients)
      names(x) <- c("Value","SE","T","P")
      x$df <- max(summary(model)$df)
      i1 <- confint(model)
      i2 <- confint(model, level = 0.9)
      x$n <- length(model$fitted.values)
    } else if(cls == "lme") {
      x <- as.data.frame(summary(model)$tTable)
      names(x) <- c("Value","SE","df", "T","P")
      i1 <- nlme::intervals(model, which = "fixed")[[1]]
      i2 <- nlme::intervals(model, which = "fixed", level = 0.9)[[1]]
      if(is.null(level)) x$n <- nrow(summary(model)$groups) else x$n <- length(unique(model$data[,level]))
    } else if(cls == "lmer") {
      if(class(model) == "glmerMod") {
        x <- as.data.frame(summary(model)$coefficients)
        names(x) <- c("Value", "SE", "T","P")
        x$note <- "Test statistic is z not T"
      } else if (grepl("Test", class(model))) {
        x <- suppressMessages(as.data.frame(lmerTest::summary(model)$coefficients))
        if(any(names(x) == "df")) {
          names(x) <- c("Value", "SE", "df","T","P")
        } else names(x) <- c("Value", "SE", "T")
      } else {
        x <- as.data.frame(summary(model)$coefficients)
        names(x) <- c("Value", "SE", "T")
      }
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
    x <- x[,c("Parameter","Value","SE","df","T","P","n","CI.95","sig.95","CI.90","sig.90","note")]
  } else if(type == "anova") {
    if(cls == "lme") {
      x <- anova(model)
      names(x) <- c("df_num","df_den","F","P")
      x <- cbind(Parameter = row.names(x), x)
    } else if(grepl("mer", class(model))) {
      x <- anova(model)
      names(x) <- c("df","SS","MS","F")
      x <- cbind(Parameter = row.names(x), x)
    }
  }
  if(!is.null(analysis)) x <- cbind(Analysis = analysis, x)
  return(x)
}
