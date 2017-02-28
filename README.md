# Belief UpdatR

[![Build Status](https://travis-ci.org/kleinschmidt/beliefupdatr.svg?branch=master)](https://travis-ci.org/kleinschmidt/beliefupdatr)
[![codecov](https://codecov.io/gh/kleinschmidt/beliefupdatr/branch/master/graph/badge.svg)](https://codecov.io/gh/kleinschmidt/beliefupdatr)


Belief updating for phonetic adaptation in R.

## Install

### Prerequisites

Parts of this package depend on having a working developer environment.  You may
already have these, but if not:

* __Mac OS__: install XCode (from the App Store)
* __Windows__: install [RTools](http://cran.r-project.org/bin/windows/Rtools/)
* __Linux__: you probably know what you're doing.

The
[RStan installation instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) 
walk through this process in more detail.

Finally, install `devtools` if it's not already installed:

```r
install.packages('devtools')
```

### From github

Second, install this package from github.

```r
devtools::install_github('kleinschmidt/beliefupdatr', args = "--preclean")
```

If you want to build the vignettes (long-form documentation), there are
additional dependencies.  To install this package with optional dependencies
(`Suggests` in the `DESCRIPTION` file) and build the vignettes:

```r
devtools::install_github('kleinschmidt/beliefupdatr', args = "--preclean",
                         dependencies=TRUE, build_vignettes=TRUE)
```

This will take a bit longer, because the vignettes run some simulations.

## Getting help

Functions are documented in the usual way:

```r
help(package='beliefupdatr')
?beliefupdatr::d_nix2_predict

# or
library(beliefupdatr)
?d_nix2_predict
```

Long-form documentation is available in the form of vignettes:

```r
vignette(package='beliefupdatr')
vignette('conjugate-updating')
```

Or to open a list of vignettes in your web browser:

```r
browseVignettes(package='beliefupdatr')
```

## Contents

* `/exec` has source for Stan models. You'll need [`rstan`](http://mc-stan.org/interfaces/rstan) to run them. Construct the full path to the model file like `file.path(devtools::inst('beliefupdatr'), 'exec', '<model.stan>')`.
* `/R/nix.R`: helper functions for the Normal-Chi^-2 distribution which is a conjugate prior on a Normal distribution with unknown mean and variance.

