library(dplyr)
library(purrr, warn.conflicts=FALSE)

context("Belief updating")

set.seed(1)

d <- data_frame(x = rnorm(100, 3, 1),
                c = sample(c('a', 'b'), 100, replace=TRUE),
                trial = 1:100)
pl0 <- list(a = nix2_params(0, 1, 10, 10), b = nix2_params(6, 1, 10, 10))

d1 <- belief_update(d, 'x', 'c', pl0)
d2 <- belief_update_batch(d, 'x', 'c', 'trial', seq(10, 100, by=10), pl0)

test_that("Updated beliefs are in column `beliefs`", {
  expect_that(d1, has_names(c(names(d), 'beliefs')))
  expect_that(d2, has_names(c('trial', 'beliefs')))
})

test_that("Beliefs have names of starting beliefs", {
  d1$beliefs %>% map(names) %>% unique() %>% unlist() %>% expect_equal(names(pl0))
  d2$beliefs %>% map(names) %>% unique() %>% unlist() %>% expect_equal(names(pl0))
})

test_that("Incremental and batch updating are equivalent", {
  with(left_join(d2, d1, by='trial'), expect_equal(beliefs.x, beliefs.y))
})

test_that("Batch updating checks column names", {
  expect_error(belief_update_batch(d, 'xx'), 'd does not have name xx')
  expect_error(belief_update_batch(d, 'x', 'cc'), 'd does not have name cc')
  expect_error(belief_update_batch(d, 'x', 'c', 't'), 'd does not have name t')
})
