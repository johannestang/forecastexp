forecastexp
=====

An R package for conducting macroeconomic forecasting experiments.

Content
-----------

This package contains a collection of code needed to replicate the
material in some of my papers. *Note*: Despite being collected in
a package the code has been written specifically for these papers and
contain very little error checking, so if you choose to use the code for
other purposes do so with caution.  

Usage
-----------
The package requires the __macrods__ package, both are easily installed using devtools:
```r
library(devtools)
install_github('johannestang/macrods')
install_github('johannestang/forecastexp')
library(forecastexp)
```

In addition the following packages should be installed : _glmnet_, _lars_, _pryr_, _quantreg_, _xtable_ and _zoo_, all available from CRAN.

The function `replicate_swjbes` gives an example of how to use the package. It
replicates part of the results in Stock, J. and M. Watson (2002): "Forecasting using 
principal components from a large number of predictors," Journal of the 
American Statistical Association, 97, 1167–1179.

For further examples see [sparsedi-rep](https://github.com/johannes/tang/sparsedi-rep) and
[l1factor-rep](https://github.com/johannes/tang/l1factor-rep).
