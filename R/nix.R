
## Functions for normal-inverse-chi-squared distribution and updating

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
  assert_that()
  assert_that(n>0 || p$nu>0 && p$kappa>0)
  within(list(), {
    nu <- p$nu + n
    kappa <- p$kappa + n
    mu <- (p$mu*p$kappa + xbar*n) / (kappa)
    sigma2 <- (p$nu*p$sigma2 + ss + n*p$kappa/(n+p$kappa) * (p$mu-xbar)^2) / nu
  })  
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
  
assertthat::onfailure(is_nix2_params) <- function(call, env) {
  paste0(deparse(call$x), ' are not valid Normal-Chi^-2 parameters')
}
