#' @import dplyr purrr tidyr
#' @importFrom reshape2 melt
NULL

check_training_data <- function(training, category, group) {
  training[[category]] <- training[[category]] %>% as.factor() %>% droplevels()
  training[[group]] <- training[[group]] %>% as.factor() %>% droplevels()
  return(training)
}

test_counts <- function(training, test, cue, category, response, group) {
  cat_levels <- training[[category]] %>% levels()
  subj_levels <- training[[group]] %>% levels()

  test[[response]] <- test[[response]] %>% factor(levels=cat_levels)
  test[[group]] <- test[[group]] %>% factor(levels=subj_levels)
  
  if (any(is.na(test[[response]]))) {
    stop("Levels mismatch between responses in test and categories in training")
  }
  if (any(is.na(test[[group]]))) {
    stop("Groups mismatch between training and test")
  }

  test_counts <- test %>%
    group_by_(group, cue, response) %>%
    tally() %>%
    spread_(key=response, value='n', fill=0) %>%
    ungroup()

  return(test_counts)
}

# Support functions for conjugate belief updating

#' Convert data to stan input for conjugate prior inference
#'
#' See exec/conj_id_lapsing_fit.stan for data format expected
#'
#' @param training Data frame with training data: 'cue', 'category', and
#'   'subject'.
#' @param test Data frame with test data: 'cue', 'response', and 'subject'
#' @param cue (Quoted) name of columns in training and test which have the cue
#'   values.
#' @param category (Quoted) name of column in training data with the correct
#'   category labels.
#' @param response (Quoted) name of column in test data with responses
#' @param group (Quoted) name of column with group identifiers
#'
#' @return a list with elements: x, y, z (training cues, subjects, and
#'   categories), n, m, l (number of observations, subject, and categories),
#'   x_test, y_test, z_test, and n_test.
#'
#' @examples
#' \dontrun{
#' supunsup::supunsup_clean %>%
#'   filter(supCond == 'unsupervised') %>%
#'   prepare_data_conj_infer_prior(training=., test=.,
#'                                 'vot', 'trueCat', 'respCat', 'subject')
#' }
#' @export
prepare_data_conj_infer_prior <- function(training, test,
                                          cue, category, response, group) {

  training <- check_training_data(training, category, group)
  test_counts <- test_counts(training, test, cue, category, response, group)

  within(list(), {
    x <- training[[cue]]
    y <- as.numeric(training[[group]])
    z <- as.numeric(training[[category]])

    n <- length(x)                        # num training observations
    m <- length(unique(z))                # num categories
    l <- length(unique(y))                # num groups

    x_test <- test_counts[[cue]]
    y_test <- as.numeric(test_counts[[group]])
    z_test_counts <-
      test_counts %>%
      select_(.dots=levels(training[[category]])) %>%
      as.matrix()

    n_test <- length(x_test)
  })
}

## generate category-by-subject matrix of sufficient stats
training_ss_matrix <- function(training, groupings, cue, ...) {
  training_ss <-
    training %>%
    group_by_(.dots=groupings) %>%
    summarise_each_(funs=funs(...), cue)

  stats <- names(list(...))

  map(stats, ~ reshape2::acast(training_ss, as.list(groupings), value.var=.x)) %>%
    set_names(stats)
}

#' @export
prepare_data_conj_suff_stats_infer_prior <- function(training, test, cue,
                                                     category, response, group) {

  training <- check_training_data(training, category, group)
  test_counts <- test_counts(training, test, cue, category, response, group)

  categories <- training[[category]] %>% levels()

  ## need category-by-subject/group matrices of sufficient stats
  training %>% 
    training_ss_matrix(list(category, group), cue,
                       xbar = mean, n = length, xsd = sd) %>%
    within({
      m <- dim(xbar)[1]
      l <- dim(xbar)[2]
      
      x_test <- test_counts[[cue]]
      y_test <- as.numeric(test_counts[[group]])
      z_test_counts <-
        test_counts %>%
        select_(.dots=levels(training[[category]])) %>%
        as.matrix()

      n_test <- length(x_test)
    })
    
  
}


# rename dimensions to make melting easier
rename_dims <- function(x, var, new_names) {
  names(dimnames(x[[var]])) <- new_names
  return(x)
}

# helper function to melt a mult-dimensional array of samples into a df
melt_samples <- function(samples, varname) {
  reshape2::melt(samples[[varname]], value.name=varname) %>%
    tbl_df
}


cond_from_idx <- function(i, lev) factor(lev[i], levels = lev)

#' Extract sampled prior beliefs from fitted model to a data_frame
#'
#' @param fitted Fitted stan model for inferring prior
#' @param cat_levels Labels for the categories.
#' @param samples Optional: Samples extracted from fitted model with
#'   \code{rstan::extract}, will override fitted if provided.
#'
#' @return a data_frame, with one row per sample and columns for
#'   \code{iterations} (sample number), \code{category}, \code{mu_0},
#'   \code{sigma_0}, \code{kappa_0}, and \code{nu_0}.
#' 
#' @export
extract_prior_samples <- function(fitted, cat_levels,
                                  samples=rstan::extract(fitted)) {

  varnames <- c('mu_0', 'sigma_0', 'kappa_0', 'nu_0')

  s <- samples %>%
    rename_dims('mu_0', c('iterations', 'cat_num')) %>%
    rename_dims('sigma_0', c('iterations', 'cat_num'))

  map(varnames, melt_samples, samples=s) %>%
    reduce(inner_join) %>%
    mutate(category = cond_from_idx(cat_num, cat_levels)) %>%
    select(-cat_num)

}

#' Extract sampled updated beliefs from fitted model to data_frame
#'
#' @inheritParams extract_prior_samples
#' @param group_levels Labels for groups/conditions (for updated)
#'
#' @return a data_frame, with one row per sample and columns for
#'   \code{iterations} (sample number), \code{category}, \code{group},
#'   \code{mu_n}, \code{sigma_n}, \code{kappa_n}, and \code{nu_n}.
#' 
#' @export
extract_updated_samples <- function(fitted, cat_levels, group_levels,
                                    samples=rstan::extract(fitted)) {

  varnames <- c('mu_n', 'sigma_n', 'kappa_n', 'nu_n')

  s <- samples %>%
    rename_dims('mu_n', c('iterations', 'cat_num', 'group_num')) %>%
    rename_dims('sigma_n', c('iterations', 'cat_num', 'group_num')) %>%
    rename_dims('kappa_n', c('iterations', 'cat_num', 'group_num')) %>%
    rename_dims('nu_n', c('iterations', 'cat_num', 'group_num'))

  map(varnames, melt_samples, samples=s) %>%
    reduce(inner_join) %>%
    mutate(category = cond_from_idx(cat_num, cat_levels),
           group = cond_from_idx(group_num, group_levels)) %>%
    select(-group_num, -cat_num)

}

#'
#'
updated_samples_predict_lhood <- function(updated_samples, x) {
  
  y <- updated_samples %>%
    by_row(~ nix2_params(.$mu_n, .$sigma_n^2, .$kappa_n, .$nu_n) %>%
             d_nix2_predict(x=x, p=.) %>%
             data_frame(x=x, pred=.)) %>%
    select(iterations, category, group, .out) %>%
    unnest(.out)

}