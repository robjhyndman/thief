## Development version
  * Fixed warning and incorrect prediction intervals in reconcilethief() when
    the forecast horizon is one period at the annual level (#19)
  * Deprecated the m argument of thief(). The seasonal period is taken from
    frequency(y).
  * Fixed several bugs: thief() with comb = "shr" or "sam" and naive or snaive
    models; reconcilethief() with forecasts not ordered finest to coarsest;
    length-1 aggregates plotted as zero; tsaggregates() with an aggregatelist
    that omits 1.

## V 0.3 (24 January 2018)
  * Improved handling of forecastfunction argument in thief()
  * Bug fixes
  * pkgdown documentation

## V 0.2 (7 September 2016)
  * First version on CRAN

## V 0.1 (22 August 2016)
  * First version on github
