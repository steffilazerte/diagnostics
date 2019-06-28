context("Transformations")

test_that("trans_plot", {
  trans_plot(m_lme)
  trans_plot(m_lmer)
})
