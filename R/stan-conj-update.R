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

  training[[category]] <- training[[category]] %>% as.factor() %>% droplevels()
  training[[group]] <- training[[group]] %>% as.factor() %>% droplevels()

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
      select_(.dots=cat_levels) %>%
      as.matrix()

    n_test <- length(x_test)
  })
}


## rename dimensions to make melting easier
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

#' Extract prior samples as a data frame
#' 
#' @export
extract_prior_samples <- function(samples, cat_levels) {

  varnames <- c('mu_0', 'sigma_0', 'kappa_0', 'nu_0')

  s <- samples %>%
    rename_dims('mu_0', c('iterations', 'cat_num')) %>%
    rename_dims('sigma_0', c('iterations', 'cat_num'))

  map(varnames, melt_samples, samples=s) %>%
    reduce(inner_join) %>%
    mutate(category = cond_from_idx(cat_num, cat_levels))

}

