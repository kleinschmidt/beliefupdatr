#' @import assertthat
#' @import purrr
#' @importFrom stats dt var
NULL

#' Conjugate updating of Normal-Chi^-2 parameters
#'
#' This uses the parametrization from
#' \href{https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf}{Murphy (2008;
#' University of British Columbia)}.
#'
#' @param x vector of observations to update with (optional if summary stats
#' are passed)
#' @param p named list of distribution parameters (mu, sigma2, nu, and kappa)
#' @param n number of observations (overrides x)
#' @param xbar observed mean (overrides x)
#' @param ss observed sum of squares (overrides x and s2)
#' @param s2 observed sample variance (ss / (n-1); overrides x)
#'
#' @return Updated parameter vector.
#' @export
nix2_update <- function(x=NA, p=list(nu=0, kappa=0, mu=0, sigma2=0), 
                        n=length(x), 
                        xbar=mean(x),
                        ss=s2*(n-1),
                        s2=var(x)) {
  assert_that(is_nix2_params(p))
  ## assert_that(! is.na(x) && is.na(ss) && is.na(s2))
  assert_that(n>0 || p$nu>0 && p$kappa>0)

  if (is.na(ss)) ss <- 0

  within(list(), {
    nu <- p$nu + n
    kappa <- p$kappa + n
    mu <- (p$mu*p$kappa + xbar*n) / (kappa)
    sigma2 <- (p$nu*p$sigma2 + ss + n*p$kappa/(n+p$kappa) * (p$mu-xbar)^2) / nu
  })
}

#' Wrap parameters into a list for updating etc.
#'
#' @param mu Mean
#' @param sigma2 Variance
#' @param kappa Mean prior pseudocount
#' @param nu Variance prior pseudocount
#'
#' @export
nix2_params <- function(mu, sigma2, kappa, nu) {
  p <- list(mu=mu, sigma2=sigma2, kappa=kappa, nu=nu)
  assert_that(is_nix2_params(p))
  return(p)
}

#' Check parameter vector for Normal-Chi^2
#'
#' Valid parameter vector has fields mu (expected mean), sigma2 (expected
#' variance), kappa (mean pseudocount), and nu (variance pseudocount/degrees of
#' freedom).
#'
#' Kappa, nu, and sigma2 must be non-negative.
#'
#' @param p Named list to check
#'
#' @export
is_nix2_params <- function(p) {
  assert_that(has_name(p, 'nu'))
  assert_that(has_name(p, 'kappa'))
  assert_that(has_name(p, 'mu'))
  assert_that(has_name(p, 'sigma2'))

  assert_that(p$nu >= 0)
  assert_that(p$kappa >= 0)
  assert_that(p$sigma2 >= 0)
  TRUE
}
  
assertthat::on_failure(is_nix2_params) <- function(call, env) {
  paste0(deparse(call$x), ' are not valid Normal-Chi^-2 parameters')
}

#' Predictive distribution for Normal-Chi^-2 distribution
#'
#' This is the probability distribution of observations drawn from a normal
#' distribution whose mean and variance follow the Normal-Chi^-2 distribution
#' with the specified parameters. It works out to be a Student's t distribution.
#'
#' @param x Point to evaluate likelihood at. Can be a vector.
#' @param p Parameters of Normal-Chi^-2 distribution.
#' @param log Whether to return log likelihood (default: FALSE)
#'
#' @export
d_nix2_predict <- function(x, p, log=FALSE) {
  assert_that(is_nix2_params(p))
  scale <- with(p, sqrt((1+kappa) * sigma2 / kappa))
  log_p <- dt((x-p$mu)/scale, df = p$nu, log=TRUE) - log(scale)
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}

#' Classification function for Normal-Chi^-2 distributions
#'
#' This computes the posterior probability of data under different Normal-Chi^-2
#' distributions.
#'
#' @param x Vector of points to classify
#' @param ps List of category Norma-Chi^-2 params
#' @param category (Optional) Which categories' posteriors to return. If
#'   missing, all are returned
#'
#' @return A length(x)-by-category matrix of posterior probabilities
#'
#' @export
nix2_classify <- function(x, ps, category=1:length(ps)) {
  map(ps, d_nix2_predict, x=x, log=FALSE) %>%
    do.call(cbind, .) %>%
    apply(1, function(x) x / sum(x)) %>%
    t() %>%
    `[`(, category)
}

## MUCH MUCH SLOWER:
nix2_classify_2 <- function(x, ps, category=1:length(ps)) {
  x %>%
    map(function(onex) map_dbl(ps, d_nix2_predict, x=onex, log=FALSE)) %>%
    map(~ .x / sum(.x)) %>%
    do.call(rbind, .)
}

## a little slower, but numerically more stable for very small lhoods.
nix2_classify_3 <- function(x, ps, category=1:length(ps)) {
  map(ps, d_nix2_predict, x=x, log=TRUE) %>%
    do.call(cbind, .) %>%
    apply(1, function(x) exp(x - daver::log_sum_exp(x))) %>%
    t() %>%
    `[`(, category)
}
