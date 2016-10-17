library(purrr)
context("Normal-Chi^-2 distributions and updating")

set.seed(1)

test_that("Invalid NIX2 parameters are rejected", {
  nix2_params(0,0,0,0)
  throws_error(nix2_params(Inf, 0, 0, 0))
  throws_error(nix2_params(0, Inf, 0, 0))
  throws_error(nix2_params(0, 0, -1, 0))
  throws_error(nix2_params(0, 0, 0, -1))
  throws_error(nix2_params(0, 0, Inf, 0))
  throws_error(nix2_params(0, 0, 0, Inf))
})

p0 <- nix2_params(mu=0, sigma2=1, kappa=1, nu=1)

test_that("Updating works for single data point", {
  p <- nix2_update(p0, 0)
  expect_false(any(is.na(p)))
  expect_equal(p, nix2_update_one(p0, 0))
})

test_that("Incremental and batch updating equivalent", {
  x <- rnorm(1000)
  expect_equal(reduce(x, nix2_update, .init = p0),
               nix2_update(p0, x))
})

test_that("Updating with zero prior pseudocount gives sample stats", {
  n <- 10
  x <- rnorm(n, mean=3, sd=3)
  p <- nix2_update(nix2_params(0, 0, 0, 0), x)
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

test_that("Samples from NIX2 are means and variances", {
  p <- nix2_params(0, 1, 1, 1)
  mv <- r_nix2(100, p)
  expect_equal(dim(mv), c(100, 2))
  expect_equal(colnames(mv), c("mean", "variance"))
})

test_that("Samples from NIX2 have correct moments", {
  p <- nix2_params(1, 3, 10, 20)
  mv <- r_nix2(1e6, p)

  theor_means <- function(p) with(p, c(mean=mu, variance=nu/(nu-2)*sigma2))
  theor_vars <- function(p) with(p, c(mean = sigma2/kappa * nu/(nu-2),
                                      variance = 2*(nu*sigma2)^2 / ((nu-2)^2 * (nu-4))))

  theoretical_moments <- function(p) rbind(theor_means(p), theor_vars(p))
  moments <- function(mv) apply(mv, 2, function(x) c(mean(x), var(x)))

  ## I just guessed at the tolerance here...
  expect_equal(moments(mv), theoretical_moments(p), tolerance=1e-3)
})

test_that("Posterior predictive density agrees with samples", {
  p <- nix2_params(1, 3, 10, 20)
  x <- seq(-10, 10)
  post_pred <- d_nix2_predict(x, p, log=TRUE)

  mv <- r_nix2(1e6, p)
  sampled <- map_dbl(x, ~ daver::log_mean_exp(dnorm(.x,
                                                    mv[ ,1],
                                                    sqrt(mv[ ,2]),
                                                    log=TRUE)))

  ## again, just totally winging it here with the tolerance
  expect_equal(post_pred, sampled, tolerance=1e-3)
})
