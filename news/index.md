# Changelog

## thief (development version)

- Fixed warning and incorrect prediction intervals in
  [`reconcilethief()`](http://robjhyndman.github.io/thief/reference/reconcilethief.md)
  when the forecast horizon is one period at the annual level
  ([\#19](https://github.com/robjhyndman/thief/issues/19))
- Deprecated the `m` argument of
  [`thief()`](http://robjhyndman.github.io/thief/reference/thief.md).
  The seasonal period is taken from `frequency(y)`.
- [`thief()`](http://robjhyndman.github.io/thief/reference/thief.md)
  with `comb = "shr"` or `"sam"` now works with naive and snaive models
- [`reconcilethief()`](http://robjhyndman.github.io/thief/reference/reconcilethief.md)
  now gives correct results when the forecasts are not ordered from
  finest to coarsest
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  of a `tsaggregates` object no longer show length-1 aggregates as zero
- [`tsaggregates()`](http://robjhyndman.github.io/thief/reference/tsaggregates.md)
  now aggregates the first series when `aggregatelist` omits 1, and
  gives an error for an invalid `m`

## thief 0.3

CRAN release: 2018-01-24

- Improved handling of `forecastfunction` argument in
  [`thief()`](http://robjhyndman.github.io/thief/reference/thief.md)
- Bug fixes
- pkgdown documentation

## thief 0.2

CRAN release: 2016-09-07

- First version on CRAN

## thief 0.1

GitHub release: 2016-08-22

- First version on github
