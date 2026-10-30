# Plot time series aggregates

Plots all temporal aggregations of a time series

## Usage

``` r
# S3 method for class 'tsaggregates'
plot(x, series = "all", ...)

# S3 method for class 'tsaggregates'
autoplot(object, series = "all", ...)
```

## Arguments

- x:

  `tsaggregates` object, produced by
  [`tsaggregates`](http://robjhyndman.github.io/thief/reference/tsaggregates.md).

- series:

  The indexes of the series to plot. By default, all series are plotted.

- ...:

  Other arguments passed to
  [`plot.ts`](https://rdrr.io/r/stats/plot.ts.html) or
  [`autoplot.ts`](https://pkg.robjhyndman.com/forecast/reference/autoplot.ts.html).

- object:

  `tsaggregates` object, produced by
  [`tsaggregates`](http://robjhyndman.github.io/thief/reference/tsaggregates.md).

## Author

Rob J Hyndman

## Examples

``` r
deathagg <- tsaggregates(USAccDeaths)
plot(deathagg, series=c(1,2,4,6))


library(ggplot2)
autoplot(deathagg)

```
