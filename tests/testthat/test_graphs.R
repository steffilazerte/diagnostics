context("Graphs")

test_that("ggResid", {
  lapply(m_list, function(x) expect_silent(ggResid(get(x))))
})

test_that("ggQQ", {
  lapply(m_list, function(x) expect_silent(ggQQ(get(x))))
})
