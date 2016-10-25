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
