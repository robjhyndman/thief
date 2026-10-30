# Temporal hierarchical forecasting

Takes a time series as input and produces forecasts using the temporal
hierarchical approach of Athanasopoulos et al (2016).

## Usage

``` r
thief(y, m = frequency(y), h = m * 2, comb = c("struc", "mse", "ols",
  "bu", "shr", "sam"), usemodel = c("ets", "arima", "theta", "naive",
  "snaive"), forecastfunction = NULL, aggregatelist = NULL, ...)
```

## Arguments

- y:

  Time series input

- m:

  Seasonal period

- h:

  Forecast horizon

- comb:

  Combination method of temporal hierarchies, taking one of the
  following values:

  "struc"

  :   Structural scaling - weights from temporal hierarchy

  "mse"

  :   Variance scaling - weights from in-sample MSE

  "ols"

  :   Unscaled OLS combination weights

  "bu"

  :   Bottom-up combination – i.e., all aggregate forecasts are ignored.

  "shr"

  :   GLS using a shrinkage (to block diagonal) estimate of residuals

  "sam"

  :   GLS using sample covariance matrix of residuals

- usemodel:

  Model used for forecasting each aggregation level:

  "ets"

  :   exponential smoothing, using the
      [`ets`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
      function.

  "arima"

  :   arima, using the
      [`auto.arima`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
      function.

  "theta"

  :   theta method, using the
      [`thetaf`](https://pkg.robjhyndman.com/forecast/reference/forecast.theta_model.html)
      function.

  "naive"

  :   random walk forecasts

  "snaive"

  :   seasonal naive forecasts, based on the last year of observed data.

- forecastfunction:

  User-defined function to be used instead of `usemodel`. The function
  must take a time series as the first argument, and the forecast
  horizon as the second argument. It must return an object of class
  `forecast`.

- aggregatelist:

  User-selected list of forecast aggregates to consider

- ...:

  Arguments to be passed to the time series modelling function (such as
  `ets` or `auto.arima`), or to `forecastfunction`.

## Value

forecast object.

## Details

This function computes the temporal aggregates of `y` using
[`tsaggregates`](http://robjhyndman.github.io/thief/reference/tsaggregates.md),
then calculates all forecasts using the model function specified by
`usemodel` or `forecastfunction`, and finally reconciles the forecasts
using
[`reconcilethief`](http://robjhyndman.github.io/thief/reference/reconcilethief.md).
The reconciled forecasts of `y` are returned.

## See also

[`reconcilethief`](http://robjhyndman.github.io/thief/reference/reconcilethief.md)

## Author

Rob J Hyndman and Nikolaos Kourentzes

## Examples

``` r
if (FALSE) { # \dontrun{

# Select ARIMA models for all series using auto.arima()
z <- thief(AEdemand[,12], usemodel='arima')
plot(z)

# Use your own function
ftbats <- function(y,h,...){forecast(tbats(y),h,...)}
z <- thief(AEdemand[,12], forecastfunction=ftbats)
plot(z)
} # }
```
