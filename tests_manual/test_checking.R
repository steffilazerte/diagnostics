library(lme4)
library(ggplot2)
data(Orthodont, package = "nlme")

m1 <- lmer(distance ~ age * Sex + (1|Subject), data = Orthodont)

vifmer(m1)
kappamer(m1)

ggplot(Orthodont, aes(x = age, y = distance, colour = Sex, group = Sex)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(Orthodont, aes(x = Sex, y = age, colour = Sex, group = Sex)) +
  geom_boxplot()
