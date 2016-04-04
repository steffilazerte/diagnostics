#' KMO for PCA
#'
#' @export
kmo = function( data ){
  ## FROM http://www.opensubscriber.com/message/r-help@stat.math.ethz.ch/7315408.html
  library(MASS)
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
#' @import lme4
#' @import Matrix
#' @export
vifmer <- function (fit) {
    ## adapted from rms::vif

    v <- vcov(fit)
    nam <- names(fixef(fit))

    ## exclude intercepts
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
#' Condition number test less than 10 or 20 (less than 30)
#' @import lme4
#' @export
kappamer <- function (fit,
                       scale = TRUE, center = FALSE,
                       add.intercept = TRUE,
                       exact = FALSE) {
    X <- getME(fit, name = "X")
    nam <- names(fixef(fit))

    ## exclude intercepts
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
sig.test <- function(model, group = "obs", alpha = 0.05){
    require(plyr)
    if(class(model) == "lm") {
        group = "obs"
        data <- model$model
        orig <- data.frame(n = 0, p = summary(model)$coefficients[,'Pr(>|t|)'])
    }
    if(class(model) == "lme") {
        data <- model$data
        orig <- data.frame(n = 0, p = summary(model)$tTable[,'p-value'])
    }
    if(grepl("mer", class(model))) {
        data <- model@frame
        c <- confint(model, method = "Wald", level = 1-alpha)
        c <- sign(c[,1]) == sign(c[,2])
        orig <- data.frame(n = 0, p = c[-grep(".sig",names(c))])
    }

    # Get original significance at alpha
    orig$var <- row.names(orig)
    orig <- dcast(orig, formula = n ~ var, value.var = "p")
    orig <- orig[,2:ncol(orig)]
    if(!grepl("mer", class(model))) orig <- orig < alpha

    # Get new sig at alpha with sequential omission
    tests <- data.frame()
    if(group == "obs") {
        for(i in 1:nrow(data)) {
            if(class(model) == "lm")  t <- data.frame(n = i, p = summary(update(model, data = data[-i,]))$coefficients[,'Pr(>|t|)'])
            if(class(model) == "lme") t <- data.frame(n = i, p = summary(update(model, data = data[-i,]))$tTable[,'p-value'])
            if(grepl("mer", class(model))) {
                c <- confint(update(model, data = data[-i,]), method = "Wald", level = 1-alpha)
                c <- sign(c[,1]) == sign(c[,2])
                t <- data.frame(n = i, p = c[-grep(".sig",names(c))])
            }
            t$var <- row.names(t)
            t <- dcast(t, formula = n ~ ..., value.var = "p")
            tests <- rbind(tests, t)
        }
    }
    if(group != "obs") {
        g <- unique(data[,group]) ## Get unique group names
        for(i in 1:length(g)) {
            if(class(model) == "lm")  t <- data.frame(n = g[i], p = summary(update(model, data = data[data[,group]!=g[i],]))$coefficients[,'Pr(>|t|)'])
            if(class(model) == "lme") t <- data.frame(n = g[i], p = summary(update(model, data = data[data[,group]!=g[i],]))$tTable[,'p-value'])
            if(grepl("mer", class(model))) {
                c <- confint(update(model, data = data[data[,group]!=g[i],]), method = "Wald", level = 1-alpha)
                c <- sign(c[,1]) == sign(c[,2])
                t <- data.frame(n = g[i], p = c[-grep(".sig",names(c))])
            }
            t$var <- row.names(t)
            t <- dcast(t, formula = n ~ ..., value.var = "p")
            tests <- rbind(tests, t)
        }
    }

    if(!grepl("mer", class(model))) tests[,2:ncol(tests)] <- tests[,2:ncol(tests)] < alpha
    output <- list(orig.sig =  orig,
                   new.sig = tests,
                   change.sig = ddply(tests, .(n), .fun = function(x, orig) x[,2:ncol(x)] != orig, orig = orig))

    return(output)
}

#' Run Diagnostics
#'
#' @export
diagnostic <- function(x, ID = "ID", alpha = 0.05, group = "obs", graphs = T, influence = T, multicol = T, power = F, group.me = NULL) {
  require(gridExtra)
  require(influence.ME)
  require(reshape2)
  if(class(x) != "lm") {
    require(lme4)
  } else require(car)

  if(class(x) == "lme") {
    require(nlme)
    require(nlmeU)
  }

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
