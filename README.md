# KW-ML for R

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/chkern/KWML.svg?branch=master)](https://travis-ci.com/chkern/KWML)
<!-- badges: end -->

R package for "[Boosted Kernel Weighting - Using Statistical Learning to Improve Inference from Nonprobability Samples](https://doi.org/10.1093/jssam/smaa028)"

Implements functions to compute pseudo-weights for nonprobability samples, including inverse propensity score weighting (`ipsw.lg()`), kernel weights based on logistic regression (`kw.lg()`), and kernel weights based on machine learning methods (`kw.mob()`, `kw.crf()`, `kw.gbm()`) 

### Installation

``` {.r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("chkern/KWML")
```

### Example

Calculate IPSW (`ipsw.lg()`) and KW-LG (`kw.lg()`) pseudo-weights with example data. `simu_dat` is a stacked data frame with a simulated probability and non-probability sample. `ipsw.lg()` and `kw.lg()` need a data frame, the name of the weight variable (`wt`, weights of 1 for non-prob, survey weights for prob sample), the name of the sample membership indicator (`trt`, 1 for non-prob, 0 for prob sample) and a formula for the propensity model as input.

``` {.r}
library(KWML)

ipsw_w <- ipsw.lg(simu_dat, "wt", "trt", 
                  "trt_f ~ x1+x2+x3+x4+x5+x6+x7")

kwlg_w <- kw.lg(simu_dat, "wt", "trt", 
                "trt_f ~ x1+x2+x3+x4+x5+x6+x7")$pswt
```

For the KW-ML functions (`kw.mob()`, `kw.crf()`, `kw.gbm()`), tuning parameter grids and covariate names for covariate balance calculation need to be specified additionally. 

``` {.r}
kwmob <- kw.mob(simu_dat, "wt", "trt", 
                "trt_f ~ x1+x2+x3+x4+x5+x6+x7 | x1+x2+x3+x4+x5+x6+x7",
                tune_maxdepth = c(2, 3), 
                covars = c("x1","x2","x3","x4","x5","x6","x7"))

kwcrf <- kw.crf(simu_dat, "wt", "trt", 
                "trt_f ~ x1+x2+x3+x4+x5+x6+x7",
                tune_mincriterion = c(0.95, 0.9), 
                covars = c("x1","x2","x3","x4","x5","x6","x7"))

kwgbm <- kw.gbm(simu_dat, "wt", "trt", 
                "trt ~ x1+x2+x3+x4+x5+x6+x7",
                tune_idepth = 1:3,
                tune_ntree = c(250, 500),
                covars = c("x1","x2","x3","x4","x5","x6","x7"))       
```

Select KW-ML pseudo-weights with best covariate balance.

``` {.r}
kwmob_w <- kwmob$pswt[, kwmob$best]
kwcrf_w <- kwcrf$pswt[, kwcrf$best]
kwgbm_w <- kwgbm$pswt[, kwgbm$best]
```

Compare weighted mean of y in prob sample and pseudo-weighted means in non-prob sample.

``` {.r}
sum((simu_dat$y[simu_dat$trt == 0]*simu_dat$wt)/sum(simu_dat$wt))
sum((simu_dat$y[simu_dat$trt == 1]*ipsw_w)/sum(ipsw_w))
sum((simu_dat$y[simu_dat$trt == 1]*kwlg_w)/sum(kwlg_w))

sum((simu_dat$y[simu_dat$trt == 1]*kwmob_w)/sum(kwmob_w))
sum((simu_dat$y[simu_dat$trt == 1]*kwcrf_w)/sum(kwcrf_w))
sum((simu_dat$y[simu_dat$trt == 1]*kwgbm_w)/sum(kwgbm_w))
```

### Citation 

``` {.r}
@Article{KernLiWang2020,
  title    = {Boosted Kernel Weighting - Using Statistical Learning to Improve Inference from Nonprobability Samples},
  author   = {Kern, C. AND Li, Y. AND Wang, L.},
  journal  = {Journal of Survey Statistics and Methodology},
  year     = {2020},
  doi      = {10.1093/jssam/smaa028},
  url      = {https://doi.org/10.1093/jssam/smaa028}
}
```
