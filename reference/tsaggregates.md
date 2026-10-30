# Non-overlapping temporal aggregation of a time series

Produces all temporal aggregations for frequencies greater than 1

## Usage

``` r
tsaggregates(y, m = frequency(y), align = c("end", "start"),
  aggregatelist = NULL)
```

## Arguments

- y:

  Univariate time series of class `ts`.

- m:

  Integer seasonal period

- align:

  Indicates how the aggregates are to be aligned: either with the
  `start` of the series or the `end` of the series. For forecasting
  purposes, it should be set to `end`.

- aggregatelist:

  User-selected list of aggregates to consider.

## Value

A list of time series. The first element is the series \`y\`, followed
by series with increasing levels of aggregation. The last element is the
"annual" series (i.e., the series aggregated over all seasons).

## See also

[`plot.tsaggregates`](http://robjhyndman.github.io/thief/reference/plot.tsaggregates.md)

## Author

Rob J Hyndman

## Examples

``` r
tsaggregates(USAccDeaths)
#> $Monthly
#>        Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec
#> 1973  9007  8106  8928  9137 10017 10826 11317 10744  9713  9938  9161  8927
#> 1974  7750  6981  8038  8422  8714  9512 10120  9823  8743  9129  8710  8680
#> 1975  8162  7306  8124  7870  9387  9556 10093  9620  8285  8466  8160  8034
#> 1976  7717  7461  7767  7925  8623  8945 10078  9179  8037  8488  7874  8647
#> 1977  7792  6957  7726  8106  8890  9299 10625  9302  8314  8850  8265  8796
#> 1978  7836  6892  7791  8192  9115  9434 10484  9827  9110  9070  8633  9240
#> 
#> $`2-Monthly`
#> Time Series:
#> Start = c(1973, 1) 
#> End = c(1978, 6) 
#> Frequency = 6 
#>  [1] 17113 18065 20843 22061 19651 18088 14731 16460 18226 19943 17872 17390
#> [13] 15468 15994 18943 19713 16751 16194 15178 15692 17568 19257 16525 16521
#> [25] 14749 15832 18189 19927 17164 17061 14728 15983 18549 20311 18180 17873
#> 
#> $Quarterly
#>       Qtr1  Qtr2  Qtr3  Qtr4
#> 1973 26041 29980 31774 28026
#> 1974 22769 26648 28686 26519
#> 1975 23592 26813 27998 24660
#> 1976 22945 25493 27294 25009
#> 1977 22475 26295 28241 25911
#> 1978 22519 26741 29421 26943
#> 
#> $`4-Monthly`
#> Time Series:
#> Start = c(1973, 1) 
#> End = c(1978, 3) 
#> Frequency = 3 
#>  [1] 35178 42904 37739 31191 38169 35262 31462 38656 32945 30870 36825 33046
#> [13] 30581 38116 34225 30711 38860 36053
#> 
#> $Biannual
#> Time Series:
#> Start = c(1973, 1) 
#> End = c(1978, 2) 
#> Frequency = 2 
#>  [1] 56021 59800 49417 55205 50405 52658 48438 52303 48770 54152 49260 56364
#> 
#> $Annual
#> Time Series:
#> Start = 1973 
#> End = 1978 
#> Frequency = 1 
#> [1] 115821 104622 103063 100741 102922 105624
#> 
#> attr(,"class")
#> [1] "tsaggregates"
```
