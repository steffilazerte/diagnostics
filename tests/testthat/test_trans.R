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
summary(m0 <- lmer(freq.top ~ spl + (1|region), data = temp))
summary(m1 <- lmer(freq.top.log ~ spl + (1|region), data = temp))

temp <- songs
songs$freq.top.log
trans(songs$freq.top.log, trans = "boxcox", boxcox = 10, centre = TRUE)

summary(m1 <- lmer(freq.top.log ~ spl + (1|region), data = temp))


test <- rnorm(10)
m <- min(test)
m.abs <- min(abs(test))

t <- "boxcox"; bx = 10
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m.abs, centre = TRUE), 9))

t <- "log"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m.abs, centre = TRUE), 9))

t <- "ln"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m.abs, centre = TRUE), 9))

t <- "inverse"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m, centre = TRUE), 9))

t <- "^2"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m.abs, centre = TRUE), 9))

t <- "^(1/2)"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m.abs, centre = TRUE), 9))

t <- "asn"; bx = NA
test2 <- (test - min(test)) / max(test - min(test))
m2 <- min(test2)
x <- trans(test2, t, bx)
all(round(test2, 9) == round(trans.back(x, t, bx, min = m2), 9))

test3 <- test2 - 0.5
m3 <- min(abs(test3))
x <- trans(test3, t, bx, centre = TRUE)
all(sort(round(test3, 9)) == sort(round(trans.back(x, t, bx, min = m3, centre = TRUE), 9)))

t <- "none"; bx = NA
x <- trans(test, t, bx)
all(round(test, 9) == round(trans.back(x, t, bx, min = m), 9))
x <- trans(test, t, bx, centre = TRUE)
all(round(test, 9) == round(trans.back(x, t, bx, min = m, centre = TRUE), 9))

