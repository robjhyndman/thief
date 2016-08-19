#thief: Temporal HIErarchical Forecasting
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/robjhyndman/thief/pulls.svg?style=flat)](https://github.com/robjhyndman/thief/pulls)

The R package *thief* provides methods and tools for generating forecasts at different temporal frequencies using a hierarchical time series approach.

Authors: Rob J Hyndman and Nikolaos Kourentzes

This package implements the methods described in 

[Athanasopoulos, G., Hyndman, R.J., Kourentzes, N., and Petropoulos, F. (2016) Forecasting with temporal hierarchies.](http://robjhyndman.com/working-papers/temporal-hierarchies/)


## Installation

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

This package is free and open source software, licensed under GPL (>= 2).
