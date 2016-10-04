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
