# Support functions for conjugate belief updating

#' Convert data to stan input for conjugate prior inference
#'
#' See exec/conj_id_lapsing_fit.stan for data format expected
#'
#' @param training Data frame with training data
#' @param test Data frame with test data
#' @param cue (Quoted) name of columns in training and test which have the cue
#'   values.
#' @param category (Quoted) name of column in training data with the correct
#'   category labels.
#' @param response (Quoted) name of column in test data with responses
#' @param subject (Quoted) name of column with subject identifiers
#'
#' @return a list with elements:

prepare_data_conj_infer_prior <- function(training, test,
                                          cue, category, response, subject) {
  test_counts <- test %>%
    group_by_(subject, cue, response) %>%
    tally() %>%
    spread_(key=response, value='n', fill=0) %>%
    ungroup()

  within(list(), {
    x <- training[[cue]]
    y <- as.numeric(training[[subject]])
    z <- as.numeric(training[[response]])

    n <- length(x)                        # num training observations
    m <- length(unique(z))                # num categories
    l <- length(unique(y))                # num subjects

    x_test <- test_counts[[cue]]
    y_test <- test_counts[[subject]]
    z_test_counts <-
      test_counts %>%
      select_(.dots=unique(test[[response]])) %>%
      as.matrix()

    n_test <- length(x_test)
  }
}
