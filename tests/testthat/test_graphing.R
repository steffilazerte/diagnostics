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
