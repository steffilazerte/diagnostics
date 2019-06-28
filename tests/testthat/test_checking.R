context("Check functions")

test_that("vif_mer calculates expected", {
  expect_equal(vif_mer(m_lme), c("wt" = 1.009952, "carb" = 1.009952), tolerance = 0.0000005)
  expect_equal(vif_mer(m_lmer), c("wt" = 1.009952, "carb" = 1.009952), tolerance = 0.0000005)
  expect_equal(vif_mer(m_glmer), c("wt" = 1.281948, "carb" = 1.281948), tolerance = 0.0000005)
  expect_equal(vif_mer(m_glmer2), c("spl" = 1.121321, "hab" = 1.110414, "regionkel" = 1.090627, "regionwil" = 1.208293 ), tolerance = 0.0000005)
})

test_that("vif() calculates vif for all model types", {
  expect_equal(car::vif(m_lm), uni_vif(m_lm))
  expect_equal(car::vif(m_glm), uni_vif(m_glm))
  expect_equal(car::vif(m_glm2), uni_vif(m_glm2))
  expect_equal(vif_mer(m_lme), uni_vif(m_lme))
  expect_equal(vif_mer(m_lmer), uni_vif(m_lmer))
  expect_equal(vif_mer(m_glmer), uni_vif(m_glmer))
  expect_equal(vif_mer(m_glmer2), uni_vif(m_glmer2))
})


test_that("uni_kappa calculates expected for all model types", {
  expect_equal(uni_kappa(m_lm), kappa(m_lm))
  #expect_equal(uni_kappa(m_glm), kappa(m_glm))
  #expect_equal(uni_kappa(m_glm2), kappa(m_lm))
  #expect_equal(uni_kappa(m_lme), kappa(m_lm))
  expect_equal(uni_kappa(m_lmer), kappa(m_lm))
  expect_equal(uni_kappa(m_glmer), kappa(m_lm))
  #expect_equal(uni_kappa(m_glmer2), kappa(m_glm2))
})

test_that("uni_kappa options", {
  expect_lt(uni_kappa(m_lm, scale = TRUE), uni_kappa(m_lm))
  expect_lt(uni_kappa(m_glm, scale = TRUE), uni_kappa(m_glm))
  #expect_lt(uni_kappa(m_lm, scale = TRUE), uni_kappa(m_lme))
  expect_lt(uni_kappa(m_lmer, scale = TRUE), uni_kappa(m_lmer))
  expect_lt(uni_kappa(m_glmer, scale = TRUE), uni_kappa(m_glmer))

  #expect_lt(uni_kappa(m_lm, intercept = FALSE), uni_kappa(m_lm))
  #expect_lt(uni_kappa(m_glm, intercept = FALSE), uni_kappa(m_glm))
  #expect_lt(uni_kappa(m_lm, intercept = FALSE), uni_kappa(m_lme))
  #expect_lt(uni_kappa(m_lmer, intercept = FALSE), uni_kappa(m_lmer))
  #expect_lt(uni_kappa(m_glmer, intercept = FALSE), uni_kappa(m_glmer))
})


test_that("overdispersion", {
  overdisp(m_glmer2)
  overdisp(m_glm2)
})

test_that("Influential variables", {
  expect_error(sig_test(m_lmer))
  expect_silent(sig_test(m_lmer2))
  expect_message(sig_test(m_lmer2, group = "cyl"))
})

test_that("multicol()", {
 expect_silent(multicol(m_lmer2))
})

test_that("Diagnostics", {
  expect_error(diagnostic(m_lmer2), NA)
  expect_error(diagnostic(m_lmer2, group = "cyl"), NA)

  diagnostic(m_lmer3, group = "site")

  expect_error(diagnostic(m_glmer2), NA)
})

temp <- pb[!pb$prob_vol & pb$ID!= "male06", ]
summary(m1 <- lmer(PC1_dist ~ cont + (1|ID), data = temp))
diagnostic(m1, group = "ID")

summary(m <- lmer(PC1_dist ~ cont + dist_start + (1|ID), data = pb))
diagnostic(m)
diagnostic(m, group = "ID")

## Can we specify different data subsets
summary(m <- lmer(PC1_dist ~ cont + dist_start + (1|ID), data = pb[!pb$prob_vol, ]))
diagnostic(m)
diagnostic(m, group = "ID")

d <- getData(m)

data(Orthodont)

songs <- read.csv("~/Projects/AAA Archive/Chorus - BCCH/Data/Datasets/bcch_pre_mean.txt")
songs$spl <- songs$spl.orig

m0 <- lm(freq ~ spl, data = songs)
diagnostic(m0)

m1 <- nlme::lme(freq ~ spl, random = ~1|region, data = songs)

m0 <- lmer(freq ~ spl + (1|region), data = songs)
m1 <- lmer(freq.log ~ spl + (1|region), data = songs)
diagnostic(m0)
diagnostic(m1)
ggQQ(m1, level = "region")

m0 <- glmer(cbind(n.top, n.bottom) ~ spl + (1|region), data = songs, family = binomial)

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

m <- lmer(distance ~ age + Sex + (1|Subject), data = Orthodont)
diagnostic(m, group = "Subject")

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

#########
## MOCH PB
#############

#' ## Load Libraries
library(seewave)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(diagnostics)

#' ## Load in Data
moch <- read.csv("~/Projects/Playback - MOCH/Data/Datasets/pb_moch_final.csv") %>%
  select(name, f.ID, f.region, f.hab, f.spl, pb.order, pb.hab.c, pb.hab, pb.spl, pb.ID, pb.region, pb.freq.l, first.rxn, start.dist, rxn.dist, rxn.song) %>%
  group_by(f.ID) %>%
  mutate(f.spl = meandB(f.spl)) %>%
  ungroup() %>%
  mutate(pb.hab.c = factor(pb.hab.c, levels = c("rural","urban"), labels = c("Rural","Urban")),
         pb.spl.c = factor(ifelse(pb.spl < median(pb.spl, na.rm = TRUE), "Quiet", "Noisy"), levels = c("Quiet", "Noisy")),
         pb.order = factor(pb.order, levels = c("U/R","R/U"), labels = c("Urban/Rural","Rural/Urban")),
         treat = factor(interaction(pb.hab.c, pb.order),
                        levels = levels(interaction(pb.hab.c, pb.order)),
                        labels = c("R.UR", "U.UR", "R.RU", "U.RU")))%>%
  filter(!(f.ID %in% c("MCC34", "MCC36")))

# /* --------------------- */
#' ## Setting up contrasts
# /* --------------------- */
contrasts(moch$treat) <- matrix(c(c(-1, 1,  0, 0),  # R vs. U in UR order
                                  c( 0, 0, -1, 1),  # R vs. U in RU order
                                  c( -1/2, -1/2, 1/2, 1/2)) # UR vs. RU
                                , ncol = 3)

#' Centre spl so that comparing to mean makes sense
moch$f.spl.orig <- moch$f.spl
moch$f.spl <- scale(moch$f.spl, scale = F)


#' ### Models
#'
round(summary(m0 <- lmer(first.rxn ~ treat + f.spl + start.dist + (1|f.ID), data = moch))$coefficients, 3)
round(summary(m1 <- lmer(first.rxn ~ treat + f.spl + (1|f.ID), data = moch))$coefficients, 3)
round(summary(m2 <- lmer(first.rxn ~ treat + (1|f.ID), data = moch))$coefficients, 3)

diagnostic(m0)
