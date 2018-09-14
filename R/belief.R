#' @import assertthat purrr
#' @importFrom dplyr %>% mutate summarise summarise_at group_by group_by_ ungroup tally select select_
#' @importFrom tidyr spread spread_ gather gather_ nest unnest
#' @importFrom utils tail
#' @importFrom rlang sym syms expr
NULL


#' Update beliefs from observations stored in a data frame
#'
#' @param d Data frame with cue values and categories
#' @param cue Quoted name of column of data frame with cue values
#' @param categories Quoted name of column with known category labels.
#' @param beliefs Named list of starting belief NIX^2 parameters. Names must
#'   match unique values of \code{d[[categories]]}
#'
#' @return The input data frame with parameter values updated after each
#'   observation in the \code{beliefs} columns
#'
#' @seealso \code{\link{belief_update_batch}} for batch belief updating, which
#'   is more efficient.
#'
#' @examples
#' d <- data.frame(x = rnorm(100, 3),
#'                 c = sample(c('a', 'b'), 100, replace=TRUE))
#' belief_update(d, 'x', 'c',
#'               list(a = nix2_params(0, 0, 10, 10),
#'                    b = nix2_params(0, 0, 10, 10)))
#'
#' @export
belief_update <- function(d, cue, categories, beliefs) {
  d$beliefs <- accumulate(transpose(d[c(cue, categories)]),
                          update_cue_category,
                          .init = beliefs) %>%
    tail(nrow(d))
  d
}

# Helper function to update beliefs for combination of cue and category
update_cue_category <- function(params, cue_category) {
  params[[ cue_category[[2]] ]] <- nix2_update_one(params[[ cue_category[[2]] ]],
                                                   cue_category[[1]])
  params
}


#' Update beliefs in a data frame in batches
#'
#' @inheritParams belief_update
#' @param trials Quoted name of column with trial numbers to define batches
#' @param at Trial numbers to get updated beliefs
#'
#' @return A data_frame with updated NIX2 parameters in \code{beliefs} for the
#'   values of \code{d[[trials]]} from \code{at}. \code{d$beliefs} is a list
#'   column where each element is a named list of category belief parameters
#'
#' @seealso \code{\link{belief_update}} for fully incremental,
#'   observation-by-observation belief updating.
#'
#' @examples
#' d <- data.frame(x = rnorm(100, 3),
#'                 c = sample(c('a', 'b'), 100, replace=TRUE),
#'                 trial = 1:100)
#' belief_update_batch(d, 'x', 'c', 'trial', seq(10, 100, by=10),
#'                     list(a = nix2_params(0, 0, 10, 10),
#'                          b = nix2_params(0, 0, 10, 10)))
#'
#' @export
belief_update_batch <- function(d, cue, categories, trials, at, beliefs) {
  assert_that(d %has_name% cue, d %has_name% categories, d %has_name% trials)
  assert_that(max(d[[trials]]) >= max(at))
  assert_that(is_empty(setdiff(d[[categories]], names(beliefs))))

  # cut, summarise by category, and update
  d %>%
    group_by(!!sym(trials) := cut_at(!!sym(trials), at)) %>%
    nest() %>%
    mutate(beliefs = map(data, ~ split(.[[cue]], .[[categories]]) %>%
                                 map(summary_stats)) %>%
             accumulate(update_summary_stats, .init=beliefs) %>%
             tail(length(data))) %>%
    select(-data)
}

summary_stats <- function(x) within(list(), {
  n <- length(x)
  xbar <- mean(x)
  ss <- var(x) * (n-1)
})

# accumulator for batch updating. params and sumstats are lists of parameters to
# update and corresponding category summary statistics respectively
update_summary_stats <- function(params, sumstats) {
  map2(params, sumstats[names(params)],
       ~ lift(nix2_update)(update_list(.y, p=.x)))
}

# label intervals with elements from 'at'
cut_at <- function(x, at) {
  codes <- .bincode(x, breaks=c(-Inf, at))
  at[codes]
}
