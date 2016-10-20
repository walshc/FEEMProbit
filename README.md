# FEEMProbit
An `R` package with a function to estimate the individual fixed effects expecation-maximization estimator in [Chen (2016)](http://blogs.bu.edu/mlchen/files/2016/03/JMP-March10th2016-version.pdf). This is work in progress and might not work as expected.

## Installation

```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("walshc/FEEMProbit")
```

## Example Usage
```r
df <- plm::pdata.frame(df, index = c("id", "time"))
FEEMProbit(y ~ x1 + x2, data = df)
```
