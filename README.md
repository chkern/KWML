# KW-ML for R

R package for "Boosted Kernel Weighting – Using Statistical Learning to Improve Inference from Nonprobability Samples"

Implements functions for inverse propensity-score weighting (IPSW; `ipsw.lg()`), kernel weights based on logistic regression (KW-LG; `kw.lg()`), and kernel weights based on machine learning methods (KW-MOB; `kw.mob()`, KW-CRF; `kw.crf()`, KW-GBM; `kw.gbm()`).

### Installation

### Example

Calculate IPSW and KW-LG pseudo-weights with example data. `simu_dat` is a stacked data frame with a simulated probability and non-probability sample. 

``` {.r}
library(KWML)

ipsw_w <- ipsw.lg(simu_dat, simu_dat$elig_wt, simu_dat$trt_n, 
                  "trt_n ~ w1 + w2 + w3 + w4 + w5 + w6 + w7")

kwlg_w <- kw.lg(simu_dat, simu_dat$elig_wt, simu_dat$trt_n, 
                "trt_n ~ w1 + w2 + w3 + w4 + w5 + w6 + w7")$pswt
```

Compare mean of y in prob sample, naive mean in non-prob sample, and pseudo-weighted means.

``` {.r}
mean(simu_dat$y[simu_dat$trt == 0])
mean(simu_dat$y[simu_dat$trt == 1])
sum((simu_dat$y[simu_dat$trt == 1]*ipsw_w)/sum(ipsw_w))
sum((simu_dat$y[simu_dat$trt == 1]*kwlg_w)/sum(kwlg_w))
```

### Citation 
