#' @importFrom dplyr inner_join
#' @importFrom abind abind
#' @importFrom reshape2 melt
NULL

#' Infer prior beliefs based on adaptation behavior
#'
#' This takes training and test data, the names of the columns, and (optionally)
#' the number of blocks to split the trials into, and stan parameters, and uses
#' stan to draw samples from the prior beliefs that match the behavior.
#'
#' @inheritParams prepare_data_incremental_suff_stats
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return A \code{stanfit} object with the fitted stan model.
#'
#' @seealso \code{\link{extract_prior_samples}} to get a data frame of samples
#'   from the prior belief parameters.
#' @export
infer_prior_beliefs <- function(df, cue, category, response, condition, ranefs,
                                n_blocks, test_df=df, ...) {

  walk(c(cue, category, response, condition, ranefs),
       ~ assert_that(has_name(df, .x)))

  if (missing(n_blocks)) {
    dat <- prepare_data_conj_suff_stats_infer_prior(df,
                                                    cue,
                                                    category,
                                                    response,
                                                    condition,
                                                    ranefs,
                                                    test_df)
    mod <- stanmodels[['conj_id_lapsing_sufficient_stats_fit']]
  } else {
    dat <- prepare_data_incremental_suff_stats(df,
                                               cue,
                                               category,
                                               response,
                                               condition,
                                               ranefs,
                                               n_blocks,
                                               test_df)
    mod <- stanmodels[['conj_id_lapsing_sufficient_stats_incremental_fit']]
  }

  fit <- rstan::sampling(mod, data=dat, ...)

}


# reduce training data by taking first ranef grouping (e.g., subject) within
# each condition, and normalize factors
prepare_training <- function(df, category, condition, ranefs) {
  df %>%
    group_by(!! sym(condition)) %>%
    filter(!!!map(syms(ranefs), ~ expr(!!.x == first(!!.x)))) %>%
    ungroup() %>%
    mutate_at(c(category, condition), . %>% as.factor() %>% droplevels())
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

## generate category-by-subject matrix of sufficient stats
training_ss_matrix <- function(training, groupings, cue, ...) {
  training_ss <-
    training %>%
    group_by_(.dots=groupings) %>%
    summarise_at(cue, funs(...))

  stats <- names(list(...))

  map(stats, ~ reshape2::acast(training_ss, as.list(groupings), value.var=.x)) %>%
    set_names(stats)
}

#' Convert data into format for bulk updating with sufficient statistics
#'
#' For Stan models that are based on sufficient statistics for training data
#' (count, mean, and sample standard deviation) instead of raw training
#' observations.
#'
#' @param df Data frame with adaptation data
#' @param cue (Quoted) name of columns in training and test which have the cue
#'   values.
#' @param category (Quoted) name of column in training data with the correct
#'   category labels.
#' @param response (Quoted) name of column in test data with responses
#' @param condition (Quoted) name of column with group identifiers
#' @param ranefs (Quoted) name of column(s) with random effect grouping
#'   variables (like subject IDs).
#' @param test_df (optional) data frame with test data (all data is used as test
#'   by default).
#'
#' @return A list of data for 'conj_id_lapsing_sufficient_stats_fit.stan'.
#'
#' @export
prepare_data_conj_suff_stats_infer_prior <- function(df, cue, category,
                                                     response, condition,
                                                     ranefs, test_df=df) {

  training <- prepare_training(df, category, condition, ranefs)

  test_counts <- test_counts(training, df, cue, category, response, condition)

  categories <- training[[category]] %>% levels()

  ## need category-by-subject/group matrices of sufficient stats
  training %>%
    training_ss_matrix(list(category, condition), cue,
                       xbar = mean, n = length, xsd = sd) %>%
    within({
      m <- dim(xbar)[1]
      l <- dim(xbar)[2]

      x_test <- test_counts[[cue]]
      y_test <- test_counts[[condition]] %>% as.numeric()
      z_test_counts <-
        test_counts %>%
        select_(.dots=levels(training[[category]])) %>%
        as.matrix()

      n_test <- length(x_test)
    })


}



#' Convert data into format for incremental updating with sufficient statistics
#'
#' Approximates true incremental updating by breaking data into blocks of equal
#' numbers of trials, and then using training data from previous n-0.5 blocks to
#' update beliefs for test data from block n.  Further, uses the _overall_
#' sufficient statistics in each group, just adjusting the number of _trials_ in
#' each block.  This is because (in the dataset this was originally developed
#' for), the trial order was different for each subject and so when fitting
#' aggregate data it makes sense to use the aggregate statistics.
#'
#' For now, the total number of trials in training and test need to be equal.
#'
#' @inheritParams prepare_data_conj_suff_stats_infer_prior
#' @param n_blocks Number of blocks to divide data into.
#'
#' @return A list of data for 'conj_id_lapsing_sufficient_stats_fit.stan'.
#'
#' @export
prepare_data_incremental_suff_stats <- function(df, cue, category, response,
                                                condition, ranefs, n_blocks,
                                                test_df = df) {

  blocks <-  df %>%
    select(trial) %>%
    mutate(block = ntile(trial, n_blocks)) %>%
    group_by(block) %>%
    summarise(max_trial = max(trial),
              min_trial = min(trial),
              mid_trial = (max_trial + min_trial)/2)

  # calculate overall summary statistics
  training <- prepare_training(df, category, condition, ranefs)

  # sufficient stats for each condition
  ss <- training_ss_matrix(training, groupings = list(category, condition), 
                           cue = cue,
                           xbar = mean, n = length, xsd = sd)

  # expand sufficient stats over blocks (just by adjusting the counts, as if
  # they've seen up to half the current block's worth of training data)
  ss_blocks <-
    blocks[['mid_trial']] %>%
                                        # fill n array, preserving names
    map(~ within(ss, n[] <- .x/dim(n)[1])) %>%
    transpose() %>%
    map(lift(abind), along=3) %>%       # make 3D array of block suff stats
    map(aperm, c(3,1,2))                # make block as first dimension

  # generate test counts broken out by block
  # the trick is to line up the group numbers. but the way they're generated for
  # the training data means that you can do group_num +
  test_counts_blocks <-
    test_df %>%
    mutate(block = cut(trial, breaks=c(min(trial)-1, blocks[['max_trial']]),
                                       labels=FALSE)) %>%
    group_by(block) %>%
    nest() %>%
    mutate(counts=map(data,
                      ~ test_counts(training, ., cue, category, response, condition))) %>%
    unnest(counts)

  within(ss_blocks, {
    k <- dim(xbar)[1]
    m <- dim(xbar)[2]
    l <- dim(xbar)[3]

    x_test <- test_counts_blocks[[cue]]
    y_test <- test_counts_blocks[[condition]] %>% as.numeric()
    block_test <- test_counts_blocks[['block']]
    z_test_counts <-
      test_counts_blocks %>%
      select_(.dots=levels(training[[category]])) %>%
      as.matrix()

    n_test <- length(x_test)
  }) %>%
  structure(blocks = blocks)

}


# rename dimensions to make melting easier
rename_dims <- function(x, var, new_names) {
  names(dimnames(x[[var]])) <- new_names
  return(x)
}


#' Melt a mult-dimensional array of samples into a df
#'
#' @param samples Named list of samples (like returned by rstan::extract)
#' @param varname Name of variable whose samples to melt
#'
#' @return A tidy dataframe of samples
#' @export
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

updated_samples_predict_lhood <- function(updated_samples, x) {

  y <- updated_samples %>%
    by_row(~ nix2_params(.$mu_n, .$sigma_n^2, .$kappa_n, .$nu_n) %>%
             d_nix2_predict(x=x, p=.) %>%
             data_frame(x=x, pred=.)) %>%
    select(iterations, category, group, .out) %>%
    unnest(.out)

}
