thief: Temporal HIErarchical Forecasting <img src="man/figures/logo.png" align="right" />
======================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/thief)](https://cran.r-project.org/package=thief)
[![Downloads](http://cranlogs.r-pkg.org/badges/thief)](https://cran.r-project.org/package=thief)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

The R package *thief* provides methods and tools for generating forecasts at different temporal frequencies using a hierarchical time series approach.

Authors: Rob J Hyndman and Nikolaos Kourentzes

This package implements the methods described in

[Athanasopoulos, G., Hyndman, R.J., Kourentzes, N., and Petropoulos, F. (2016) Forecasting with temporal hierarchies.](http://robjhyndman.com/publications/temporal-hierarchies/)


## Installation
You can install the **stable** version on
[R CRAN](https://cran.r-project.org/package=thief).

```s
install.packages('thief', dependencies = TRUE)
```

You can install the **development** version from
[Github](https://github.com/robjhyndman/thief)

```s
# install.packages("devtools")
devtools::install_github("robjhyndman/thief")
```

## Usage

```s
library(thief)
thief(USAccDeaths)
```

## License

This package is free and open source software, licensed under GPL 3.
