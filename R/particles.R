#' @import purrr
NULL

# Fearnhead (2004) algorithm:
#
# Each particle represents state of cluster assignments at time n =
# z_{1..n}^i.  with weight w_i.
#
# Observation at time n+1 comes in. Generate new enhanced set of M proposed
# particles. (z_{1..n}^i, j) for each possible assignment j. Updated weight is
# proportional to
#
#   w_i * p((z_{1:n}^i,j) | y_{1:n}) / p(z_{1:n}^i | y_{1:n}).
#
# Then resample N particles (original population size) from M. Two methods:
#
# Chen & Liu: For each particle at time n, draw assignment of y_{n+1},
#   (z^i_{1:n}, j) with j proportional to p( (z_{1:n}, j) | y_{1:n+1} ). Then
#   weight these new particles with
#
#   w^i_{n+1} = w^i_n p( (z_{1:n}^i, j) | y_{1:n+1} ) / p( z_n^i | y_{1:n} )
#
#   When weights are degenerate, rejuvinate by resampling with probability
#   proportional to w^i.
#
# Fearnhead & Clifford: Of teh M putative particles, Keep everything with weight
#   > threshold (1/c), and resample from remaining particles, assigning weight
#   1/c to them.
#
# Either way, each particle needs to be able to compute p(z | y), which you can
# do based on sufficient statistics for each category. This means that they need
# to hold onto the sufficient statistics of each category given z, and the
# posterior probability p(z | y), in addition to the weight (it's not enough to
# just hold onto the weight because the weights will change during resampling.
#
# The steps are
#
# 1. For each particle, generate all possible assignments of y_{n+1}.
# 2. For each of these descendents, calculate posterior p( (z_n^i,j) | y_{1:n+1} )
#     * update summary stats of cateogry j
#     * compute p( y_{1:n+1} | (z_n^i, j) ) (marginal likelihood of each category,
#       see Murphy eq 171). This can be done with nix2_update() with summary stats,
#       wrapped by a nix2_marginal_likelihood() function

#' Initialize one particle based on NIX2 parameters
#'
#' @param params NIX2 parameters
#' @param n Population size
#'
#' @return A (population of) particle(s)
#'
#' @export
init_particle <- function(params) {
  walk(params, ~ assert_that(is_nix2_params(.)))

  within(list(), {
    w <- 1
    z <- vector()
    params <- map(params, ~ list(prior = .x,
                                 sumstats = list(n=0, xbar=0, ss=0),
                                 marg_log_lhood = 0))
  })
}

#' @describeIn init_particle Initialize a population of particles
#'
#' @export
init_particles <- function(params, n) {
  replicate(n, init_particle(params), simplify=FALSE)
}

summary_stats_update <- function(ss, x) {
  within(ss, {
    n <- n+1
    delta <- x - xbar
    xbar <- xbar + delta/n
    ss <- ss + delta * (x-xbar)
    delta <- NULL
  })
}

params_update <- function(params, x) {
  within(params, {
    sumstats <- summary_stats_update(sumstats, x)
    ## use a "lifted" nix2_marginal_lhood so we can pass summary stats as list
    ## instead of unwrapping sumstats$n, $xbar, $ss.
    marg_log_lhood <- lift_dl(nix2_marginal_lhood)(p=prior, sumstats, log=TRUE)
  })
}

# create a new particle with x assigned to z and update weight
assign_data <- function(particle, x, new_z) {
  within(particle, {
    z <- c(z, new_z)
    last_ll <- params[[new_z]]$marg_log_lhood
    params[[new_z]] <- params_update(params[[new_z]], x)
    ## all the change in the log-posterior is from the updated category's
    ## likelihood (at least assuming a flat prior like we do here)
    w <- w * exp(params[[new_z]]$marg_log_lhood - last_ll)
  })
}

update_particle_chen_liu <- function(particle, x) {
  map(names(particle$params),
      ~ assign_data(particle, x, .)) %>%
    sample(1, prob=map_dbl(., "w")) %>%
    unlist(recursive=FALSE)
}

# we need to re-normalize weights because we're not computing exactly p(z|y), 
# but rather the un-normalized posterior p(y|z)p(z).  In the ratio of these, the
# normalizing constants p(y) don't cancel because the set of data points is
# different.  But the error is the same across particles so it's just a matter
# of scaling the weights to have unit sum.
normalize_weights <- function(particles) {
  sum_w <- sum(map_dbl(particles, 'w'))
  map(particles, update_list, w = ~ w / sum_w)
}

filter_chen_liu <- function(particles, x) {
  map(particles, update_particle_chen_liu, x) %>%
    normalize_weights()
}

################################################################################
# Doing something useful with particles:

#' Convert prior + sufficient stats to NIX2 params
particle_params_to_nix2 <- function(params) {
  assert_that(has_name(params, 'prior'), has_name(params, 'sumstats'))
  lift(nix2_update)(p=params$prior, params$sumstats)
}

#' Posterior predictive density functions
#'
#' Posterior predictive distribution for a single particle
#'
#' @param x Points to calculate posterior predictive density at
#' @param particle The particle
#' @param particles List of particles in population (weights must sum to 1)
#' @param log =FALSE
#'
#' @return A named list of predictive (log-)densities for each category (as a
#'   vector, one per x), or a data_frame with columns x, category, pred
#'   ([log-]density), particle_id, and weight
#'
#' @export
d_predict_particle <- function(x, particle, log=FALSE) {
  assert_that(has_name(particle, 'params'))
  map(particle$params,
      . %>%
        particle_params_to_nix2() %>%
        d_nix2_predict(p=., x=x, log=log))
}

#' @describeIn d_predict_particle Marginal posterior predictive distribution for
#'   population of particles (the weighted average of the individual particles'
#'   posterior predictive distributions)
#'
#' @export
d_predict_particles_marginal <- function(x, particles, log=FALSE) {
  ws <- map_dbl(particles, 'w')
  assert_that(all.equal(sum(ws), 1))

  preds <- map(particles, d_predict_particle, x=x, log=FALSE) %>%
    transpose() %>%                     # from list of particles to list of
                                        # categories
    map(~ lift(cbind)(.) %*% ws) %>%    # take weighted average of each
                                        # particles' prediction
    map(as.vector)

  if (log) {
    map(preds, base::log)
  } else {
    preds
  }
}

#' @describeIn d_predict_particle Posterior predictive densities for all
#'   particles in a population, organized into a data_frame.
#'
#' @export
d_predict_particles_df <- function(x, particles, log=FALSE) {

  map(particles, d_predict_particle, x=x, log=log) %>%
    map(~ as_data_frame(.) %>%
          mutate(x=x) %>%
          gather('category', 'pred', -x)) %>%
    map2(seq_along(.), ~ mutate(.x, particle_id=.y)) %>%
    map2(map_dbl(particles, 'w'), ~ mutate(.x, weight=.y)) %>%
    bind_rows()

}

#' @describeIn d_predict_particle Marginal posterior predictive distribution
#'   function for a whole population, as a data_frame
#'
#' @export
d_predict_particles_marginal_df <- function(x, particles, log=FALSE) {

  d_predict_particles_marginal(x, particles, log) %>%
    as_data_frame() %>%
    mutate(x=x) %>%
    gather('category', 'pred', -x)

}
