#' KMO for PCA
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
  } else {
    g <- lapply(unique(data[, group]), FUN = function(x) which(data[, group] == x))
  }

  # Get new sig at alpha with sequential omission
  tests <- data.frame()
  for(i in 1:length(g)) {
    if(verbose) print(paste0(i, " / ", length(g)))

    t <- cbind(n = i,
               get.table(my.update(model, data = data[-g[[i]], ]))[, c("Parameter", "Value", "P")]
    )
    tests <- rbind(tests, t)
  }

  output <- merge(tests, orig[, c("Parameter", "Value", "P")], by = "Parameter", suffixes = c("", ".orig"))
  output$diff5 <- (output$P < 0.05) != (output$P.orig < 0.05)
  output$diff10 <- (output$P < 0.10) != (output$P.orig < 0.10)
  output <- output[output$diff5 | output$diff10, ]

  return(output)
}

#' Run Diagnostics
#'
#' @import gridExtra
#' @import influence.ME
#' @import reshape2
#' @import lme4
#' @import car

#' @export
diagnostic <- function(x, ID = "ID", alpha = 0.05, group = "obs", graphs = T, influence = T, multicol = T, power = F, group.me = NULL) {

  require(nlme)
  require(nlmeU)
  ## Get residual plots and normality plots for all levels of the random variables (if present)
  if(graphs == T){
    g <- list()
    g[[1]] <- ggResid(x) + labs(title = "Residual Plot") ## Residual
    g[[2]] <- ggQQ(x) + labs(title = "QQ Normality Plot") ## General Normality
    if(class(x) == "lme" | grepl("mer", class(x))) {
      for(a in 3:(length(ranef(x))+2)) {
        if(length(ranef(x)[[a-2]][[1]]) > 2 & mean(abs(ranef(x)[[a-2]][[1]])) > 0){
          g[[length(g)+1]] <- ggQQ(x, level = a-2) + labs(title = paste0("QQ Normality Plot Level ", a-2)) ## Normality of random
        }
      }
    }
    do.call("grid.arrange", c(g, nrow = 1))  ## Plot
  }

  ## Get power
  if(power == T & class(x) == "lme") p <- Pwr(x)

  m <- x
  if(class(x) == "lme") m <- lmer(formula(paste0(deparse(x$call$fixed), " + (", gsub("\\~(+)", "\\1",deparse(x$call$random)),")")), data = x$data, control=lmerControl(optimizer="bobyqa"))

  ## Get influential observations
  if(influence == T){
    n <- "Not calculable"
    v <- "Not applicable"
    k <- "Not applicable"

    if(!any(grep("~ 1", getCall(x)))) {
      i <- sigtest(x, group = group, alpha = alpha)[[3]]
      n <- vector()
      ## Get cooks
      cooks <- vector()
      if(class(x) != "lm"){
        cks <- round(cooks.distance(influence(m, obs = ifelse(is.null(group.me), T, F), group = group.me)),3)
      } else cks <- round(influence.measures(m)$infmat[,'cook.d'],3)

      for(a in 1:nrow(i)) {
        if(any(i[a,-1] == T)) {
          n <- c(n,paste0(i$n[a], " (",paste0(names(i[,-1])[which(i[a,-1]==T)], collapse = ", "),")"))
          cooks <- c(cooks, cks[a])
        }
      }

      if(length(n) == 0) n <- "None"
    } else cooks <- "N/A"
  }

  ## Get Multiple colinearity
  if(multicol == T & !any(grep("~ 1", getCall(x)))){
    if(class(x) != "lm") {
      v <- vif.mer(m)
      k <- kappa.mer(m)
    } else {
      if(length(m$coefficients) > 2) v <- vif(m) else v <- NA
      if(!is.null(dim(v))) v <- v[,1]
      k <- kappa(m)
    }
  } else {v <- NA; k <- NA}

  ## Display
  if(influence == T) if(length(n) > 0 & class(x) == "lm") print("NOTE: Influence identifies rows in a data frame from which all missing variables are removed.")
  if(influence == T) print(paste0("Influence: ", paste0(n, collapse = "; ")))
  if(influence == T) print(paste0("Cooks Distances: ", paste0(cooks, collapse = "; ")))
  if(multicol == T) {
    print(paste0("VIF: ", paste(names(v), round(v, digits = 1), sep = " ")))
    print(paste0("Kappa: ", round(k,digits = 1)))
  }
  if(power == T & class(x) == "lme") {print(paste0("Power alpha = 0.05")); print(p)}

}

#' @export
get.table <- function(m, analysis = NULL, type = "summary", level = NULL, pca = 1){
  if(type == "summary") {
    if(class(m) == "lm"){
      x <- as.data.frame(summary(m)$coefficients)
      names(x) <- c("Value","SE","T","P")
      x$df <- max(summary(m)$df)
      i1 <- confint(m)
      i2 <- confint(m, level = 0.9)
      x$n <- length(m$fitted.values)
    } else if(class(m) == "lme") {
      x <- as.data.frame(summary(m)$tTable)

      names(x) <- c("Value","SE","df", "T","P")
      i1 <- intervals(m, which = "fixed")[[1]]
      i2 <- intervals(m, which = "fixed", level = 0.9)[[1]]
      if(is.null(level)) x$n <- nrow(summary(m)$groups) else x$n <- length(unique(m$data[,level]))
    } else if(grepl("mer", class(m))) {
      x <- as.data.frame(summary(m)$coefficients)
      if(class(m) == "glmerMod") {
        names(x) <- c("Value", "SE", "T","P")
        x$note <- "Test statistic is z not T"
      } else if (grepl("Test", class(m))) {
        x <- as.data.frame(lmerTest::summary(m)$coefficients)
        names(x) <- c("Value", "SE", "df","T","P")
      } else {
        names(x) <- c("Value", "SE", "T")
      }
      i1 <- confint(m, method = "Wald")
      i1 <- i1[-grep(".sig",row.names(i1)),]
      i2 <- confint(m, method = "Wald", level = 0.9)
      i2 <- i2[-grep(".sig",row.names(i2)),]
      if(is.null(level)) x$n <- nrow(m@frame) else x$n <- length(unique(m@frame[,level]))
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
    if(class(m) == "lme") {
      x <- anova(m)
      names(x) <- c("df_num","df_den","F","P")
      x <- cbind(Parameter = row.names(x), x)
    } else if(grepl("mer", class(m))) {
      x <- anova(m)
      names(x) <- c("df","SS","MS","F")
      x <- cbind(Parameter = row.names(x), x)
    }
  } else if(type == "pca" & class(m) == "prcomp") {
    x <- round(as.data.frame(m$rotation)[,1:pca],3)
    for(i in 1:pca) x[,paste0("Var",i)] <- round(summary(m)$importance[2,i],3)
    x <- cbind(Parameter = row.names(x),x)
  }
  if(!is.null(analysis)) x <- cbind(Analysis = analysis, x)
  return(x)
}
