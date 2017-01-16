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

# Fearnhead (2004) resampling scheme:
# 
# 1. Generate all M putative particles from assignments of data to category for
#    each particle and update weights.
# 
# 2. Based on new weights, find cut-off $1/c$ such that
#    $\sum_{i=1}^M min(cw_i, 1) = N$.
#
# 3. Keep all the particles with weights > 1/c, and resample from the remaining
#    particles. (using stratified sampling so that each is sampled at most once)

# update one particle with all possible assignments
update_particle_fearnhead <- function(particle, x) {
  map(names(particle$params), ~ assign_data(particle, x, .))
}

putative_particles <- function(particles, x) {
  particles %>%
    map(update_particle_fearnhead, x) %>%
    unlist(recursive=FALSE) %>%
    normalize_weights()
}


# algorithm is from Clifford, Carpenter, and Fearnhead (1999), as described in
# Fearnhead and Clifford (2003, Appendix B)
#
# I THINK THIS IS WRONG. (but probably okay for these purposes where N << M)
#
# does not always return teh requested number of samples (fewer sometimes, never
# more).  But this only seems to be a big problem when you've got size ≈
# length(x), which isn't going to happen often (I think).
sample_stratified <- function(x, size, prob) {
  K <- sum(prob) / size
  U <- runif(1, 0, K)
  include <- rep(FALSE, length(x))
  for (i in seq_along(x)) {
    U <- U - prob[i]
    if (U < 0) {
      include[i] <- TRUE
      U <- U + K
    }
  }
  x[include]
}

# resample particles so that there's at most N, according to the algorithm from
# Fearnhead and Clifford (2003).  This preserves all the particles that have
# weight above some cutoff, and resamples the rest and sets the weights to
# 1/cutoff
resample_fearnhead <- function(particles, N) {
  # sort by decreasing weight
  M <- length(particles)
  if (M <= N) return(particles)

  # should be 1 but just in case...
  tot <- particles %>% map_dbl("w") %>% sum()

  # find cutoff of weights to include.  return value is c, such that
  #   sum(min(c*w, 1)) = N.
  #
  # strategy is to go through weights from smallest to largest, to find the first
  #   w[i] where sum(w[1:i]) / w[i] < (N-(M-i)).  Then 1/c is between w[i-1] and
  #   w[i], and we can easily solve for it with sum(w[1:i-1]) * c + M - (i-1) = N.
  #   if you run out of w's before you find the right i, then c = 1/N.
  #
  # actually: more efficient to go largest to smallest. then once we find first i
  #   where sum(w[1:i])/w[i] >= (N-(M-i)), solve for sum(w[1:i])*c + M-i = N.
  #
  # (that's more efficient because it's more likely that there's none that are kept
  particles <- sort_by(particles, ~ -.x[['w']])
  # (i will be the first weight that is _not_ automatically carried over.)
  for (i in seq_along(particles)) {
    # so 1:(i-1) contriubtes 1 each, and the rest contribute
    # less than sum(w[i:end]) / w[i] (1/c will be <= w[i])
    if (tot / particles[[i]][['w']] + i-1 >= N) {
      # solve for c in sum(w[i:end]) * c + i-1 = N
      cutoff <- (N - (i-1)) / tot
      break
    }
    # incrementally update sum:
    tot <- tot - particles[[i]][['w']]
  }

  # not sure if this can happen, but make sure we return at most N
  if (i > N) return(particles[1:N])
  
  # keep particles[i:i-1], and resample N-(i-1) from particles[i:end]
  new_particles <- 
    particles[i:M] %>%
    sample_stratified(., N-(i-1), map_dbl(., 'w')) %>%
    map(update_list, w = 1/cutoff) %>%
    c(particles[1:(i-1)], .)

  return(new_particles)
}

filter_fearnhead <- function(particles, x, N=length(particles)) {
  particles %>%
    map(update_particle_fearnhead, x) %>%
    unlist(recursive=FALSE) %>%
    normalize_weights() %>%
    resample_fearnhead(N)
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
