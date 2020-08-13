# KW-ML for R

R package for "Boosted Kernel Weighting – Using Statistical Learning to Improve Inference from Nonprobability Samples"

Implements functions for inverse propensity-score weighting (`ipsw.lg()`), kernel weights based on logistic regression (`kw.lg()`), and kernel weights based on machine learning methods (`kw.mob()`, `kw.crf()`, `kw.gbm()`).

### Installation

### Example

Calculate IPSW (`ipsw.lg()`) and KW-LG (`kw.lg()`) pseudo-weights with example data. `simu_dat` is a stacked data frame with a simulated probability and non-probability sample. Both functions need the data frame, the name of the weight variable (`wt`), the name of the sample membership indicator (`trt`) and a formula for the propensity model as input.

``` {.r}
library(KWML)

ipsw <- ipsw.lg(simu_dat, "wt", "trt", 
                "trt_f ~ x1 + x2 + x3 + x4 + x5 + x6 + x7")

kwlg <- kw.lg(simu_dat, "wt", "trt", 
              "trt_f ~ x1 + x2 + x3 + x4 + x5 + x6 + x7")$pswt
```

For the KW-ML functions (`kw.mob()`, `kw.crf()`, `kw.gbm()`), tuning parameter grids and covariate names for covariate balance calculation need to be specified additionally. 

``` {.r}
kwmob <- kw.mob(simu_dat, "wt", "trt", 
               "trt_f ~ x1+x2+x3+x4+x5+x6+x7 | x1+x2+x3+x4+x5+x6+x7",
               c(2, 3), 
               c("x1","x2","x3","x4","x5","x6","x7"))

kwcrf <- kw.crf(simu_dat, "wt", "trt", 
               "trt_f ~ x1+x2+x3+x4+x5+x6+x7",
               c(0.95, 0.9), 
               c("x1","x2","x3","x4","x5","x6","x7"))
               
kwgbm <- kw.gbm(simu_dat, "wt", "trt", 
               "trt ~ x1+x2+x3+x4+x5+x6+x7",
               1:3,
               c(250, 500),
               c("x1","x2","x3","x4","x5","x6","x7"))               
```

Compare weighted mean of y in prob sample and pseudo-weighted means in non-prob sample.

``` {.r}
sum((simu_dat$y[simu_dat$trt == 0]*simu_dat$wt)/sum(simu_dat$wt))
sum((simu_dat$y[simu_dat$trt == 1]*ipsw)/sum(ipsw))
sum((simu_dat$y[simu_dat$trt == 1]*kwlg)/sum(kwlg))
```

### Citation 
