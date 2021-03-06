#' @import assertthat
#' @importFrom stats dt var
NULL

#' Conjugate updating of Normal-Chi^-2 parameters
#'
#' This uses the parametrization from
#' \href{https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf}{Murphy (2008;
#' University of British Columbia)}.
#'
#' @param p named list of distribution parameters (mu, sigma2, nu, and kappa)
#' @param x vector of observations to update with (optional if summary stats
#' are passed)
#' @param n number of observations (overrides x)
#' @param xbar observed mean (overrides x)
#' @param ss observed sum of squares (overrides x and s2)
#' @param s2 observed sample variance (ss / (n-1); overrides x)
#'
#' @return Updated parameters
#' @export
nix2_update <- function(p=list(nu=0, kappa=0, mu=0, sigma2=0),
                        x=NA,
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
    mu <- (p$mu*p$kappa + xbar*n) / kappa
    sigma2 <- (p$nu*p$sigma2 + ss + n*p$kappa/kappa * (p$mu-xbar)^2) / nu
  })
}


#' Conjugate updating of Normal-Chi^-2 parameters for a single observation
#'
#' This uses the parametrization from
#' \href{https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf}{Murphy (2008;
#' University of British Columbia)}.
#'
#' This is specialized for a single observation, and about 25% faster than
#' \code{nix2_update}. Use, for instance, in \code{reduce} or \code{accumulate}.
#'
#' @param p named list of distribution parameters (mu, sigma2, nu, and kappa)
#' @param x vector of observations to update with (optional if summary stats
#' are passed)
#'
#' @return Update parameters
#' @export
nix2_update_one <- function(p=list(nu=0, kappa=0, mu=0, sigma2=0), x) {
  assert_that(is_nix2_params(p))
  within(list(), {
    nu <- p$nu + 1
    kappa <- p$kappa + 1
    mu <- (p$mu * p$kappa + x) / kappa
    sigma2 <- (p$nu * p$sigma2 + p$kappa/kappa * (p$mu - x)^2) / nu
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

  assert_that(all(p$nu >= 0))
  assert_that(all(p$kappa >= 0))
  assert_that(all(p$sigma2 >= 0))
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
    lift(cbind)() %>%
    apply(1, function(x) x / sum(x)) %>%
    t() %>%
    `[`(, category)
}

#' Marginal likelihood function
#'
#' Compute the marginal likelihood of data given NIX^2 prior, or \eqn{p(x) =
#' \int \mathrm d \sigma_2 \mathrm d\mu p(x | \mu, \sigma^2) p(\mu, \sigma^2)}.
#'
#' Based on equation 171 from
#' \href{https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf}{Murphy (2008;
#' University of British Columbia)}.
#'
#' @inheritParams nix2_update
#' @param log Default FALSE.
#'
#' @return The overal marginal likelihood of all the observations x (or summary
#'   statistics if provided).
#' @export
nix2_marginal_lhood <- function(p,
                                x=NA,
                                n=length(x),
                                xbar=mean(x),
                                ss=s2*(n-1),
                                s2=var(x),
                                log = FALSE) {
  pn <- nix2_update(p, n=n, xbar=xbar, ss=ss)
  L <- nix2_log_norm_const(pn) - nix2_log_norm_const(p) - (n/2)*logpi
  ifelse(log, L, exp(L))
}

logpi <- log(pi)

#' Normalizing constant for NIX^2 distribution
#'
#' The terms of \eqn{p(\mu, \sigma^2 | \mu_0, \sigma^2_0, \nu, \kappa)} that
#' don't involve \eqn{\mu} or \eqn{\sigma^2}.  Useful in calculating the
#' marginal likelihood of data, where the mean and variance are integrated out.
#'
#' @param p NIX^2 parameters
nix2_log_norm_const <- function(p) {
  with(p, lgamma(nu/2) - 0.5 * (log(kappa) + nu * (log(nu) + log(sigma2))))
}

#' Sample from NIX2 distribution
#'
#' This first draws the variance from a (scaled) inverse chi squared
#' distribution, and then the mean from a normal distribution.
#'
#' @param n Number of samples to draw
#' @param p NIX2 parameters
#'
#' @return An n x 2 matrix with sampled means in the first column and variances
#'   in the second
#' @export
r_nix2 <- function(n, p) {
  assert_that(is_nix2_params(p))

  vars <- p$nu * p$sigma2 / rchisq(n, df = p$nu)
  means <- rnorm(n, mean = p$mu, sd = sqrt(vars/p$kappa))

  cbind(mean=means, variance=vars)
}

#' Scaled inverse chi-squared density
#'
#' We compute this based on the chi-squared density function.  Our random
#' variable y is obtained from a chi-squared(nu) variable x via
#'
#'   y = g(x) = nu * sigma^2 * 1/x
#'   x = g^-1(y) = nu*sigma^2 * 1/y
#'
#' Thus, the density function f_Y(y) is obtained from f_X = dchisq via
#'
#'   f_Y(y) = d/dy g^-1(y) * f_X(g^-1(y))
#'          = y^-2 * nu*sigma^2 * dchisq(nu*sigma^2/y, nu)
#'
#' @param y Point to calculate density for
#' @param sigma2 Scale parameter
#' @param nu Degrees of freedom
#' @param log Default: FALSE
#'
#' @return The density of y under the scaled-inverse-chi-squared distribution
#' @export
d_scaled_inv_chisq <- function(y, sigma2, nu, log=FALSE) {
  log_d <- -2*log(y) + log(nu) + log(sigma2) + dchisq(nu*sigma2/y, df=nu, log=TRUE)
  if (log) log_d
  else exp(log_d)
}

#' Probability density function for NIX2 distribution
#'
#' @param mv Matrix of mean-variance pairs, means in the first column and
#'   variances in the second.
#' @param p NIX2 parameters.
#' @param log FALSE to return probability (default), TRUE for log-probability.
#'
#' @return The (log-)probability density for each row of \code{mv}.
#' @export
d_nix2 <- function(mv, p, log=FALSE) {
  assert_that(is_nix2_params(p))
  if (is.null(dim(mv))) {
    assert_that(length(mv) == 2)
    mv <- matrix(mv, ncol=2)
  } else {
    assert_that(dim(mv)[2] == 2)
  }

  log_d <-
    with(p,
         dnorm(mv[, 1], mean=mu, sd=sqrt(mv[, 2] / kappa), log=TRUE) +
           d_scaled_inv_chisq(mv[, 2], sigma2=sigma2, nu=nu, log=TRUE))

  if (log) log_d
  else exp(log_d)

}
