library(dplyr, warn.conflicts=FALSE)


context("Prepare stan data")

d <- supunsup::supunsup_clean %>% filter(supCond == 'unsupervised')

test_that('Data converts for incremental update', {

  d_inc <- prepare_data_incremental_suff_stats(training = d %>%
                                                 group_by(bvotCond) %>%
                                                 filter(subject == first(subject)),
                                               test = d,
                                               cue = 'vot',
                                               category = 'trueCat',
                                               response = 'respCat',
                                               group = 'bvotCond',
                                               n_blocks = 10)

  expect_that(d_inc, has_names(c('xbar', 'n', 'xsd', 'l', 'm', 'k',
                                 'n_test', 'z_test_counts', 'block_test',
                                 'y_test', 'x_test'), ignore.order=TRUE))

  d_inc[c('xbar', 'n', 'xsd')] %>%
    map(dim) %>%
    walk(expect_equal, with(d_inc, c(k, m, l)))

  expect_that(attributes(d_inc), has_names(c('names', 'blocks')))
})
