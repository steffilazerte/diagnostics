
test.plot <- function(model, level = "all", type = "QQ") {
  par(mfrow=c(1, 2), mar = c(2,2,2,2))
  if(type == "QQ") {
    if(level == "all"){
      qqnorm(residuals(model), main = "QQ Plot")
      qqline(residuals(model))
    } else {
      qqnorm(ranef(model)[[level]][[1]], main = "QQ Plot")
      qqline(ranef(model)[[level]][[1]])
    }
  } else if (type == "R") plot(fitted(model), residuals(model))
  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
  vp1 <- plotViewport(c(0,0,0,0))
  if(type == "QQ") g <- ggQQ(model, level = level) else if(type == "R") g <- ggResid(model)
  print(g, vp = vp1)
}

library(grid)
library(gridBase)
library(lme4)
library(nlme)

data(Orthodont)

m.lm1 <- lm(distance ~ age, data = Orthodont)

m.lme1 <- lme(distance ~ age, random = ~1|Sex, data = Orthodont)
m.lme2 <- lme(distance ~ age, random = ~1|Sex/Subject, data = Orthodont)

m.lmer1 <- lmer(distance ~ age + (1|Sex), data = Orthodont)
m.lmer2 <- lmer(distance ~ age + (1|Sex/Subject), data = Orthodont)
m.lmer3 <- lmer(distance ~ age + (1|Sex) + (1|Subject), data = Orthodont)

m.glmer1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)

#' ## Test function errors
ggQQ(m.lm1)
ggResid(m.lm1)

#' ### One random factor
ggResid(m.lme1)
ggQQ(m.lme1)
ggQQ(m.lme1, level = "Sex")

ggResid(m.lmer1)
ggQQ(m.lmer1)
ggQQ(m.lmer1, level = "Sex")

#' ### Nested random factors
ggResid(m.lme2)
ggQQ(m.lme2)
ggQQ(m.lme2, level = "Subject")

ggResid(m.lmer2)
ggQQ(m.lmer2)
ggQQ(m.lmer2, level = "Subject")

#' ### Crossed random factors
ggResid(m.lmer3)
ggQQ(m.lmer3)
ggQQ(m.lmer3, level = "Sex")
ggQQ(m.lmer3, level = "Subject")

#' ## GLME
ggResid(m.glmer1)
ggQQ(m.glmer1)
ggQQ(m.glmer1, level = "herd")

#' ## Test Plotting
test.plot(m.lm1, type = "R")
test.plot(m.lm1)

test.plot(m.lme1, type = "R")
test.plot(m.lme1)
test.plot(m.lme1, level = "Sex")

test.plot(m.lme2, type = "R")
test.plot(m.lme2)
test.plot(m.lme2, level = "Subject:Sex")

test.plot(m.lmer1, type = "R")
test.plot(m.lmer1)
test.plot(m.lmer1, level = "Sex")

test.plot(m.lmer2, type = "R")
test.plot(m.lmer2)
test.plot(m.lmer2, level = "Subject:Sex")

test.plot(m.lmer3, type = "R")
test.plot(m.lmer3)
test.plot(m.lmer3, level = "Sex")
test.plot(m.lmer3, level = "Subject")


# Test return W
ggQQ(m.lm1) # W = 0.986
ggQQ(m.lmer1) # W = 0.988
ggQQ(m.lmer2) # W = 0.926

ggQQ(m, plot = FALSE) # W = 0.986
ggQQ(m.lmer1, plot = FALSE) # W = 0.988
ggQQ(m.lmer2, plot = FALSE) # W = 0.926
