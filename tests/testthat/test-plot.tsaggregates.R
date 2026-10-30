test_that("plot.tsaggregates draws without error", {
  agg <- tsaggregates(USAccDeaths)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(agg))
  expect_no_error(plot(agg, series = c(1, 2, 4, 6)))
  expect_no_error(plot(agg, series = 3, main = "Quarterly"))
})

test_that("plot.tsaggregates validates series", {
  agg <- tsaggregates(USAccDeaths)
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(agg, series = 7), "series should be between 1 and 6")
  expect_error(plot(agg, series = 0), "series should be between 1 and 6")
  big <- tsaggregates(ts(rnorm(2 * 720), frequency = 720))
  expect_error(plot(big), "Too many series")
})

test_that("autoplot.tsaggregates returns a ggplot", {
  skip_if_not_installed("ggplot2")
  agg <- tsaggregates(USAccDeaths)
  p <- ggplot2::autoplot(agg)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_s3_class(ggplot2::autoplot(agg, series = c(1, 6)), "ggplot")
  expect_error(ggplot2::autoplot(agg, series = 10), "series should be between")
})

test_that("mts.from.tsaggregates picks axis labels from the frequency", {
  agg <- tsaggregates(USAccDeaths)
  expect_equal(mts.from.tsaggregates(agg, "all")$args$xlab, "Year")
  expect_equal(mts.from.tsaggregates(agg, "all", xlab = "X")$args$xlab, "X")
  daily <- tsaggregates(ts(1:28, frequency = 7))
  expect_equal(mts.from.tsaggregates(daily, "all")$args$xlab, "Week")
})
