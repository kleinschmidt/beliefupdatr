library(purrr)
context("Normal-Chi^-2 distributions and updating")

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
  p <- nix2_update(0, p0)
  expect_false(any(is.na(p)))
})

test_that("Incremental and batch updating equivalent", {
  x <- rnorm(1000)
  expect_equal(reduce(x, ~ nix2_update(x=.y, p=.x), .init = p0),
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


test_that("Classification always between 0 and 1", {

  ps <- list(nix2_params(0, 1, 10, 10),
             nix2_params(3, 1, 10, 10))
  x <- seq(-4, 7, by=.1)
  clas <- nix2_classify(x, ps)

  expect_true(all(clas >= 0))
  expect_true(all(clas <= 1))
})

test_that("Classification of each point adds up to 1", {
  ps <- list(nix2_params(0, 1, 10, 10),
             nix2_params(3, 1, 10, 10))
  x <- seq(-4, 7, by=.1)
  clas <- nix2_classify(x, ps)
  
  expect_equal(apply(clas, 1, sum), rep(1, dim(clas)[1]))
})

test_that("Classification function is ~logistic with high confidence", {
  mus <- c(0,1)
  ps <- map(mus, nix2_params, sigma2=1, kappa=1e9, nu=1e9)
  x <- seq(-1, 2, by=0.1)
  clas <- nix2_classify(x, ps, category=1)

  # slope is -(mu[1]-mu[2])/sigma^2, intercept is (mu[1]^2-mu[2]^2)/(2*sigma^2)
  # (cf. Feldman et al., 2009, p. 781)
  clas_pred <- (1 + exp(x-0.5))^-1
  expect_equal(clas, clas_pred)
})
          
