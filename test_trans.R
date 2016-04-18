library(diagnostics)
library(lme4)

data(Orthodont, package = "nlme")


m1 <- lmer(distance ~ age * Sex + (1|Subject), data = Orthodont)

trans.plot(m1, find.best = TRUE)

trans.plot(m1)
trans.plot(m1, t = "boxcox")


songs <- read.csv("~/Projects/Chorus - BCCH/Data/Datasets/bcch_exp_mean_difference.txt")
head(songs$freq.top)
summary(m0 <- lmer(freq.top ~ spl + (1|region), data = songs))

temp <- songs
temp$freq.top <- trans(temp$freq.top, trans = "boxcox", boxcox = 2)
head(temp$freq.top)
summary(update(m0, data = temp))

trans.plot(m0)
trans.plot(m0, verbose = TRUE)

trans.plot(m0, t = "boxcox", boxcox = seq(1, 5, 0.5))
trans.plot(m0, t = "boxcox", boxcox = seq(1, 5, 0.5), verbose = TRUE)

trans.plot(m0, find.best = TRUE)
trans.plot(m0, find.best = TRUE, verbose = TRUE)

head(songs$freq.top)
head(trans(songs$freq.top, trans = "log"))
head(trans(songs$freq.top, trans = "inverse"))
head(trans(songs$freq.top, trans = "^3"))

temp$freq.top.log <- trans(temp$freq.top.log, trans = "boxcox", boxcox = 10)



summary(m1.a <- lmer(freq.top.log ~ spl + (1|region), data = temp))

temp <- songs
summary(m0.b <- lmer(freq.top ~ spl + (1|region), data = temp))
summary(m1.b <- lmer(freq.top.log ~ spl + (1|region), data = temp))
