# Belief UpdatR

Belief updating for phonetic adaptation in R.

This is really quite disorganized at this point, use at your own risk.

## Install

```r
devtools::install_github('kleinschmidt/beliefupdatr')
```

## Contents

* `/exec` has source for Stan models. You'll need [`rstan`](http://mc-stan.org/interfaces/rstan) to run them. Construct the full path to the model file like `file.path(devtools::inst('beliefupdatr'), 'exec', '<model.stan>')`.
* `/R/nix.R`: helper functions for the Normal-Chi^-2 distribution which is a conjugate prior on a Normal distribution with unknown mean and variance.
