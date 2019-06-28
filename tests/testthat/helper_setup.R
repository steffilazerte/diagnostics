
library(magrittr)
songs <- tibble::as_tibble(read.csv("~/Projects/AAA Archive/Chorus - MOCH/Docs/Manuscripts/Data/data_pre_chorus.csv"))
songs <- dplyr::mutate(songs, ident = 1:length(p.s))

m_list <- c("m_lm", "m_glm", "m_glm2", "m_lme", "m_lmer", "m_lmer2", "m_glmer", "m_glmer2")

m_lm <- lm(mpg ~ wt + carb, data = mtcars)
m_glm <- glm(vs ~ wt + carb, data = mtcars, family = binomial)
m_glm2 <- glm(cbind(p.s, p.c) ~ spl + hab + region, family = binomial, data = songs)
m_lme <- nlme::lme(mpg ~ wt + carb, random = ~1|cyl, data = mtcars)
m_lmer <- lme4::lmer(mpg ~ wt + carb + (1|cyl), data = mtcars)
m_lmer2 <- lmerTest::lmer(mpg ~ wt + carb + (1|cyl), data = mtcars)
m_lmer3 <- lmerTest::lmer(p.t ~ spl + hab + (1|site), data = songs)
m_glmer <- lme4::glmer(vs ~ wt + carb + (1|cyl), data = mtcars, family = binomial)
m_glmer2 <- lme4::glmer(cbind(p.s, p.c) ~ spl + hab + region + (1|ident), family = binomial, data = songs)
