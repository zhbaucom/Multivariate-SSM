# Multivariate-SSM
Functions to create and estimate multivariate state space models

## The Model

In the project “SSM for Repeated Measures”, we established that the
Bayesian LLT was superior to commonly used LMEMs in fitting longitudinal
data. This project extends that method to allow for multivariate
Bayesian LLT, meaning, multiple outcomes. The important gains from this
approach is to increase power, estimate inter-relatedness in tests after
accounting for predictors of interest, and lastly being able to compare
effects across different outcomes.

To measure outcome inter-relatedness we allow for a subject's observation
errors and latent cognition processes to be correlated. Estimating these
correlation process can be used to detect larger underlying cognition
constructs. To compare linear effects, we propose an in-Gibb’s sampler
standardization to ensure outcomes are on the same scale. This process
relies on accurate estimation of observation error correlation. One cannot simply
standardize outcomes beforehand because it will not take into account
the variance coming from the latent cognition process.

## Fitting the Model

Once again, we utilize the Bayesian Gibb’s sample as it has been shown
superior to the other likelihood based counterparts.

## How to Use

``` r
source("functions/JBLLT.R")
JBLLT(
  method = "os", data, outcomes, predictors, timevar, id, 
  initialization = "Bayes", numits = 1000, burn = 500, 
  silence = FALSE, seed = NULL, 
  numitInit = 1500, burnInit = 500, seed.init = NULL
)
```

Here we provide all data for fitting. The `outcomes`, `predictors`,
`timevar`, and `id` should all be characters corresponding to columns in
the `data`. Notice that you must clarify if their is allowed correlation in the observation equation `method = "o"`, in the state equation `method = "s"`, or both `method = "os"` (default).
