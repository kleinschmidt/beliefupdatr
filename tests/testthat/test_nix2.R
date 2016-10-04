context("Normal-Chi^-2 distributions and updating")


p0 <- nix2_params(mu=0, sigma2=1, kappa=1, nu=1)

test_that("Updating works for single data point", {
  p <- nix2_update(0, p0)
  expect_false(any(is.na(p)))
})

test_that("Incremental and batch updating equivalent", {
  x <- rnorm(1000)
  expect_equal(purrr::reduce(x, ~ nix2_update(x=.y, p=.x), .init = p0),
               nix2_update(x, p0))
})

test_that("Updating with zero prior pseudocount gives sample stats", {
  n <- 10
  x <- rnorm(n, mean=3, sd=3)
  p <- nix2_update(x, nix2_params(0, 0, 0, 0))
  expect_equal(p$mu, mean(x))
  expect_equal(p$sigma2, var(x)*(n-1)/n)
  expect_equal(p$kappa, n)
  expect_equal(p$nu, n)
})

test_that("Posterior predictive is nearly normal for high confidence", {
  x <- seq(-3, 3, by=.1)
  p <- nix2_params(0, 1, 1e9, 1e9)
  expect_equal(dnorm(x), d_nix2_predict(x, p))
})

