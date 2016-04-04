/*
 * Fit multinomial response data using a belief-updating model to infer prior
 * parameters.  A normal-inverse-chi-squared prior is used, and it's assumed
 * that the subject knows the true labels of all the input stimuli (e.g.,
 * doesn't model any uncertainty in classification.
 *
 * This version has a lapse rate parameter (probability of random guessing)
 *
 * Dave Kleinschmidt
 * Apr 2015
 */

data {
  int n;                        // number of training observations
  real x[n];                    // training observations
  int m;                        // number of categories
  int z[n];                     // categories of each training observation
  int l;                        // number of subjects
  int y[n];                     // subject labels
  int n_test;                   // number of test trials
  real x_test[n_test];          // locations of test trials
  int y_test[n_test];           // subject labels for test trials
  int z_test_counts[n_test,m];  // responses for test trials
}

transformed data {
  real xbar[m,l];               // mean
  real ss[m,l];                 // sum of squares
  real n_cat[m,l];              // number of obs per category

  xbar <- rep_array(0, m,l);
  ss <- rep_array(0, m,l);
  n_cat <- rep_array(0, m,l);

  // running mean/sum-of-squares calculations
  for (j in 1:n) {
    real delta;
    int cat;
    int subj;
    cat <- z[j];
    subj <- y[j];
    n_cat[cat,subj] <- n_cat[cat,subj] + 1;
    delta <- x[j] - xbar[cat,subj];
    xbar[cat,subj] <- xbar[cat,subj] + delta / n_cat[cat,subj];
    ss[cat,subj] <- ss[cat,subj] + delta * (x[j] - xbar[cat,subj]);
  }

}

parameters {
  // these are all shared across subjects (same prior beliefs):
  real<lower=0> kappa_0;        // prior pseudocount for mean
  real<lower=0> nu_0;           // prior pseudocount for sd
  real mu_0[m];                 // prior expected mean
  real<lower=0> sigma_0[m];     // prior expected standard deviation
  real<lower=0, upper=1> lapse_rate;
}

transformed parameters {
  // updated beliefs depend on input/subject
  real mu_n[m,l];                 // updated expected mean
  real<lower=0> kappa_n[m,l];     // updated mean pseudocount
  real<lower=0> sigma_n[m,l];     // updated expected sd
  real<lower=0> nu_n[m,l];        // updated sd pseudocount
  real<lower=0> t_scale[m,l];     // scale parameter of predictive t distribution
  simplex[m] p_test_conj[n_test];
  vector[m] log_p_test_conj[n_test];

  // update NIX2 parameters according to conjuate updating rules are taken from
  // Murphy (2007)
  for (cat in 1:m) {
    for (subj in 1:l) {
      kappa_n[cat,subj] <- kappa_0 + n_cat[cat,subj];
      nu_n[cat,subj] <- nu_0 + n_cat[cat,subj];
      mu_n[cat,subj] <- (mu_0[cat] * kappa_0 + xbar[cat,subj] * n_cat[cat,subj]) / kappa_n[cat,subj];
      sigma_n[cat,subj] <- sqrt((nu_0*sigma_0[cat]^2 +
                                 ss[cat,subj] +
                                 (n_cat[cat,subj]*kappa_0)/(kappa_n[cat,subj]) *
                                   (mu_0[cat] - xbar[cat,subj])^2
                                 ) /
                                nu_n[cat,subj]);
      t_scale[cat,subj] <- sigma_n[cat,subj] * sqrt((kappa_n[cat,subj] + 1) / kappa_n[cat,subj]);
    }
  }

  // compute category probabilities for each of the test stimuli
  for (j in 1:n_test) {
    int subj;
    subj <- y_test[j];
    // calculate un-normalized log prob for each category
    for (cat in 1:m) {
      log_p_test_conj[j,cat] <- student_t_log(x_test[j],
                                              nu_n[cat,subj],
                                              mu_n[cat,subj],
                                              t_scale[cat,subj]);
    }
    // normalize and store actual probs in simplex
    p_test_conj[j] <- exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  real n_each;
  vector[m] lapsing_probs;
  
  n_each <- n / (m*l);
  lapsing_probs <- rep_vector(lapse_rate / m, m);
  
  // need to calculate category probabilities for each test trial
  kappa_0 ~ normal(0, n_each*4);
  nu_0 ~ normal(0, n_each*4);

  mu_0 ~ normal(0, 100);
  sigma_0 ~ uniform(0, 100);

  for (i in 1:n_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }
  
}

generated quantities {
  // simplex[m] p_test_semi_conj[n_test];
  // simplex[m] p_test_sampled[n_test];
  // real ll;

  // ll <- 0;
  // for (i in 1:n_test) {
  //   ll <- ll + multinomial_log(z_test_counts[i], p_test_conj[i]);
  // }
  
}
