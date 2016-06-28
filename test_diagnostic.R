library(lme4)
library(nlme)
library(lmerTest)
library(diagnostics)

data(Orthodont)

songs <- read.csv("~/Projects/Chorus - BCCH/Data/Datasets/bcch_pre_mean.txt")
songs$spl <- songs$spl.orig

m0 <- lm(freq ~ spl, data = songs)
diagnostic(m0)

m0 <- lmer(freq ~ spl + (1|region), data = songs)
m1 <- lmer(freq.log ~ spl + (1|region), data = songs)
diagnostic(m0)
diagnostic(m1)
ggQQ(m1, level = "region")

diagnostic(m0, influence = FALSE)
diagnostic(m0, influence = FALSE, multicol = FALSE)

songs <- read.csv("~/Projects/Chorus - BCCH/Data/Datasets/bcch_exp_mean_difference.txt")
summary(m0 <- lmer(freq.change.log ~ spl + (1|region), data = songs))
diagnostic(m0)

songs$freq.change[1:5]
trans(songs$freq.change, "boxcox", 3, centre = TRUE)[1:5]

m.lm1 <- lm(distance ~ age, data = Orthodont)
m.lm2 <- lm(distance ~ age + Subject, data = Orthodont)

m.lme1 <- lme(distance ~ age, random = ~1|Sex, data = Orthodont)
m.lme2 <- lme(distance ~ age, random = ~1|Sex/Subject, data = Orthodont)

m.lmer1 <- lmer(distance ~ age + (1|Sex), data = Orthodont)
m.lmer2 <- lmer(distance ~ age + (1|Sex/Subject), data = Orthodont)
m.lmer3 <- lmer(distance ~ age + (1|Sex) + (1|Subject), data = Orthodont, REML = TRUE)

m.lmer1.2 <- lmer(distance ~ age + Subject + (1|Sex), data = Orthodont)

temp <- Orthodont
temp$distance[1] <- 500
temp$age[5] <- NA
m.lmer1.3 <- lmer(distance ~ age+ (1|Sex), data = temp)

m.glmer1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)

summary(m.lmer1)
summary(as(m.lmer1, "lmerMod"))
summary(as(m.lmer1, "lmerTest"))
summary(m <- update(m.lmer1, data = Orthodont))
summary(m <- my.update(m.lmer1, data = Orthodont))
summary(as(m, "merModLmerTest"))

names(summary(my.update(m.lmer1, data = Orthodont)))
summary(as(my.update(m.lmer1, data = Orthodont), "merModLmerTest"))

sig.test(m.lm1)
sig.test(m.lme1, all = TRUE, group = "Subject")
sig.test(m.lme2)
sig.test(m.lme2, verbose = TRUE)
sig.test(m.lmer1, verbose = TRUE)
sig.test(m.lmer2, verbose = TRUE)
sig.test(m.lmer3, verbose = TRUE)
sig.test(m.glmer1, verbose = TRUE)

multicol(m.lm1)
multicol(m.lme1)
multicol(m.lme2)
multicol(m.lmer1)
multicol(m.lmer2)
multicol(m.lmer3)
multicol(m.lmer1.2)

diagnostic(m.lm1)
diagnostic(m.lme1)
diagnostic(m.lme2)
diagnostic(m.lmer1)
diagnostic(m.lmer1.2, verbose = TRUE)
diagnostic(m.lmer2)
diagnostic(m.lmer3)
