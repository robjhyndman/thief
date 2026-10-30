# thief (development version)

- Fixed warning and incorrect prediction intervals in `reconcilethief()` when the forecast horizon is one period at the annual level (#19)
- Deprecated the `m` argument of `thief()`.
  The seasonal period is taken from `frequency(y)`.
- `thief()` with `comb = "shr"` or `"sam"` now works with naive and snaive models
- `reconcilethief()` now gives correct results when the forecasts are not ordered from finest to coarsest
- `plot()` and `autoplot()` of a `tsaggregates` object no longer show length-1 aggregates as zero
- `tsaggregates()` now aggregates the first series when `aggregatelist` omits 1, and gives an error for an invalid `m`

# thief 0.3

- Improved handling of `forecastfunction` argument in `thief()`
- Bug fixes
- pkgdown documentation

# thief 0.2

- First version on CRAN

# thief 0.1

<p class="text-muted">GitHub release: 2016-08-22</p>

- First version on github
