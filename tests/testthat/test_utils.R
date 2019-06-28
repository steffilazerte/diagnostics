context("Utility functions")

test_that("is_CLASS returns correct assesment", {
  expect_true(is_lm(m_lm))
  expect_false(is_lme(m_lm))
  expect_false(is_mer(m_lm))

  expect_true(is_lm(m_glm))
  expect_false(is_lme(m_glm))
  expect_false(is_mer(m_glm))

  expect_false(is_lm(m_lme))
  expect_true(is_lme(m_lme))
  expect_false(is_mer(m_lme))

  expect_false(is_lm(m_lmer))
  expect_false(is_lme(m_lmer))
  expect_true(is_mer(m_lmer))

  expect_false(is_lm(m_glmer))
  expect_false(is_lme(m_glmer))
  expect_true(is_mer(m_glmer))
})

test_that("get_fixed() returns fixed effects", {
  expect_equal(get_fixed(m_lm), summary(m_lm)$coeff[, 1])
  expect_equal(get_fixed(m_glm), summary(m_glm)$coeff[, 1])
  expect_equal(get_fixed(m_lme), summary(m_lme)$coeff$fixed)
  expect_equal(get_fixed(m_lmer), summary(m_lmer)$coeff[, 1])
  expect_equal(get_fixed(m_glmer), summary(m_glmer)$coeff[, 1])
})

test_that("get_data() returns full model dataset", {
  expect_equal(get_data(m_lm), mtcars)
  expect_equal(get_data(m_glm), mtcars)
  expect_equal(get_data(m_lme), mtcars)
  expect_equal(get_data(m_lmer), mtcars)
  expect_equal(get_data(m_glmer), mtcars)
})

test_that("get_data() returns subset model dataset", {
  m_lm <- lm(mpg ~ wt + carb, data = mtcars[mtcars$mpg < 20, ])
  m_glm <- glm(vs ~ wt + carb, data = mtcars[mtcars$mpg < 20, ], family = binomial)
  m_lme <- nlme::lme(mpg ~ wt + carb, random = ~1|cyl, data = mtcars[mtcars$mpg < 20, ])
  m_lmer <- lme4::lmer(mpg ~ wt + carb + (1|cyl), data = mtcars[mtcars$mpg < 20, ])
  m_glmer <- lme4::glmer(vs ~ wt + carb + (1|cyl), data = mtcars[mtcars$mpg < 20, ], family = binomial)

  expect_equal(get_data(m_lm), mtcars[mtcars$mpg < 20, ])
  expect_equal(get_data(m_glm), mtcars[mtcars$mpg < 20, ])
  expect_equal(get_data(m_lme), mtcars[mtcars$mpg < 20, ])
  expect_equal(get_data(m_lmer), mtcars[mtcars$mpg < 20, ])
  expect_equal(get_data(m_glmer), mtcars[mtcars$mpg < 20, ])
})

test_that("model_update() updates model", {

  expect_silent(model_update(m_lm, data = mtcars[mtcars$mpg < 20, ]))
  expect_silent(model_update(m_glm, data = mtcars[mtcars$mpg < 20, ]))
  expect_silent(model_update(m_lme, data = mtcars[mtcars$mpg < 20, ]))
  expect_silent(model_update(m_lmer, data = mtcars[mtcars$mpg < 20, ]))
  expect_silent(model_update(m_glmer, data = mtcars[mtcars$mpg < 20, ]))
  #summary(model_update(m_lm, data = mtcars[mtcars$mpg < 20, ]))

  expect_equivalent(class(m_lm), class(model_update(m_lm, data = mtcars)))
  expect_equivalent(class(m_glm), class(model_update(m_glm, data = mtcars)))
  expect_equivalent(class(m_lme), class(model_update(m_lme, data = mtcars)))
  expect_equivalent(class(m_lmer), class(model_update(m_lmer, data = mtcars)))
  expect_equivalent(class(m_lmer2), class(model_update(m_lmer2, data = mtcars)))
  expect_equivalent(class(m_glmer), class(model_update(m_glmer, data = mtcars)))
  expect_equivalent(class(m_glmer2), class(model_update(m_glmer2, data = songs)))

  expect_equivalent(broom::tidy(m_lm), broom::tidy(model_update(m_lm, data = mtcars)))
  expect_equivalent(broom::tidy(m_glm), broom::tidy(model_update(m_glm, data = mtcars)))
  expect_equivalent(broom::tidy(m_lme), broom::tidy(model_update(m_lme, data = mtcars)))
  expect_equivalent(broom::tidy(m_lmer), broom::tidy(model_update(m_lmer, data = mtcars)))
  expect_equivalent(broom::tidy(m_lmer2), broom::tidy(model_update(m_lmer2, data = mtcars)))
  expect_equivalent(broom::tidy(m_glmer), broom::tidy(model_update(m_glmer, data = mtcars)))
  expect_equivalent(broom::tidy(m_glmer2), broom::tidy(model_update(m_glmer2, data = songs)))
  })
