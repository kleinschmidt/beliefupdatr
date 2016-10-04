context("Classifier given NIX^2 beliefs")

mus <- c(a=0,b=1)
ps <- map(mus, nix2_params, sigma2=1, kappa=1e9, nu=1e9)
x <- seq(-1, 2, by=0.1)
clas <- nix2_classify(x, ps)

test_that("Classification always between 0 and 1", {
  expect_true(all(clas >= 0))
  expect_true(all(clas <= 1))
})

test_that("Classification of each point adds up to 1", {
  expect_equal(apply(clas, 1, sum), rep(1, dim(clas)[1]))
})

test_that("Classification function is ~logistic with high confidence", {
  # slope is -(mu[1]-mu[2])/sigma^2, intercept is (mu[1]^2-mu[2]^2)/(2*sigma^2)
  # (cf. Feldman et al., 2009, p. 781)
  clas_pred <- (1 + exp(x-0.5))^-1
  expect_equal(clas[,1], clas_pred)
})

test_that("Index classification with numeric index", {
  expect_equal(nix2_classify(x, ps, category=2), clas[, 2])
  throws_error(nix2_classify(x, ps, category=-1))
})

test_that("Category names are preserved as column names in classification", {
  expect_equal(colnames(clas), names(mus))
})

test_that("Index classification with category name", {
  expect_equal(nix2_classify(x, ps, category='b'), clas[, 2])
  throws_error(nix2_classify(x, ps, category='c'))
})
