# KW-ML for R

R package for "Boosted Kernel Weighting – Using Statistical Learning to Improve Inference from Nonprobability Samples"

Implements functions for inverse propensity-score weighting (IPSW; `ipsw.lg()`), kernel weights based on logistic regression (KW-LG; `kw.lg()`), and kernel weights based on machine learning methods (KW-MOB; `kw.mob()`, KW-CRF; `kw.crf()`, KW-GBM; `kw.gbm()`).

### Installation

### Example

Calculate IPSW and KW-LG pseudo-weights with example data. `simu_dat` is a stacked data frame with a simulated probability and non-probability sample. 

``` {.r}
library(KWML)

ipsw <- ipsw.lg(simu_dat, "wt", "trt", 
                "trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7")

kwlg <- kw.lg(simu_dat, "wt", "trt", 
              "trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7")$pswt
```

Compare weighted mean of y in prob sample and pseudo-weighted means in non-prob sample.

``` {.r}
sum((simu_dat$y[simu_dat$trt == 0]*simu_dat$wt)/sum(simu_dat$wt))
sum((simu_dat$y[simu_dat$trt == 1]*ipsw)/sum(ipsw))
sum((simu_dat$y[simu_dat$trt == 1]*kwlg)/sum(kwlg))
```

### Citation 
