
# take a data frame with observations, and accumulate through nix2_update_one.

d <- 
p0 <- nix2_params(0, 1, kappa=10, nu=10)

d$p <- purrr::accumulate(d$x, nix2_update_one, .init=p0) %>% tail(nrow(d))

data_frame(x = rnorm(100, 3),
           trial = 1:length(x)) %>%
  mutate(p = purrr::accumulate(x, nix2_update_one, .init=p0) %>% tail(length(x)))



d2 <- data_frame(x = rnorm(100, 3),
                 c = sample(c("a", "b"), size=100, replace=TRUE),
                 trial = 1:length(x))

d2 %<>%
  group_by(c) %>%
  mutate(p = purrr::accumulate(x, nix2_update_one, .init=p0) %>% tail(length(x)))



#' Update beliefs in a data frame
#'
#' @param d Data frame with cue values and (optionally) categories
#' @param cue Quoted name of column of data frame with cue values
#' @param categories (Optional) name of column with known category labels.
#' @param params Starting NIX^2 parameter values.
#'
#' @return The input data frame with parameter values updated after each
#'   observation in the \code{params} columns
belief_update <- function(d, cues, categories=NULL, params) {
  q <- lazyeval::interp(~ purrr::accumulate(x, nix2_update_one, .init=params) %>%
                          tail(length(x)),
                        x = as.name(cues))
  d %>%
    group_by_(.dots=categories) %>%
    mutate_(params = q)
}

belief_update(d, 'x', NULL, p0)

# the problem with this approach
# 1. awkward to restrict which params to update
# 2. no flexibility in how to handle missing categories.
#
# better design seems to be to write a separate reducing function that takes
# named list of params, cues, and categories, and updates the matching params.
# (but how to do this efficiently isn't clear. need to zip obs and categories...)
# 


# Helper function to update beliefs for combination of cue and category
update_cue_category <- function(params, cue_category) {
  params[[ cue_category[[2]] ]] <- nix2_update_one(params[[ cue_category[[2]] ]],
                                                   cue_category[[1]])
  params
}

belief_update_2 <- function(d, cues, categories, params) {
  d$params <- purrr::accumulate(purrr::transpose(d[c(cues, categories)]),
                                update_cue_category,
                                .init = params) %>%
    tail(nrow(d))
  d
}

cc <- transpose(d[c('x', 'c')])
  
pl0 <- list('a'=p0, 'b'=p0)
d3 <- belief_update_2(d, 'x', 'c', pl0)

# surprisingly, version with transpose isn't much slower (~1ms out of 20ms)
microbenchmark(belief_update_2(d, 'x', 'c', pl0),
               belief_update(d, 'x', 'c', p0))
