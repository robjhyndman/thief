test_that("tsaggregates returns named aggregates that conserve totals", {
  agg <- tsaggregates(USAccDeaths)
  expect_s3_class(agg, "tsaggregates")
  expect_named(agg, c("Monthly", "2-Monthly", "Quarterly", "4-Monthly", "Biannual", "Annual"))
  expect_equal(agg$Monthly, USAccDeaths)
  expect_equal(frequency(agg$Annual), 1)
  expect_length(agg$Annual, 6)
  expect_equal(sum(agg$Annual), sum(USAccDeaths))
  expect_equal(frequency(agg$Quarterly), 4)
  expect_equal(sum(agg$Quarterly), sum(USAccDeaths))
})

test_that("tsaggregates names quarterly aggregates", {
  agg <- tsaggregates(ts(1:16, frequency = 4))
  expect_named(agg, c("Quarterly", "Biannual", "Annual"))
})

test_that("tsaggregates align argument controls which observations are dropped", {
  y <- ts(1:14, frequency = 4, start = 2000)
  end <- tsaggregates(y, align = "end")
  start <- tsaggregates(y, align = "start")
  # 14 observations: 3 full years; end alignment drops the first 2, start the last 2
  expect_equal(as.numeric(end$Annual), c(sum(3:6), sum(7:10), sum(11:14)))
  expect_equal(as.numeric(start$Annual), c(sum(1:4), sum(5:8), sum(9:12)))
})

test_that("tsaggregates respects aggregatelist", {
  agg <- tsaggregates(USAccDeaths, aggregatelist = c(1, 3, 12))
  expect_named(agg, c("Monthly", "Quarterly", "Annual"))
})

test_that("tsaggregates works with a non-ts input and explicit m", {
  agg <- tsaggregates(1:12, m = 4)
  expect_equal(frequency(agg[[length(agg)]]), 1)
  expect_equal(sum(agg[[length(agg)]]), sum(1:12))
})

test_that("tsaggregates errors on bad input", {
  expect_error(tsaggregates(USAccDeaths, aggregatelist = 5), "No valid factors")
  expect_error(tsaggregates(1:5, m = 0), "too short")
})

test_that("tsaggregates aggregates the first series when aggregatelist omits 1", {
  agg <- tsaggregates(USAccDeaths, aggregatelist = c(3, 12))
  expect_named(agg, c("Quarterly", "Annual"))
  expect_equal(frequency(agg$Quarterly), 4)
  expect_length(agg$Quarterly, 24)
  expect_equal(sum(agg$Quarterly), sum(USAccDeaths))
  expect_equal(frequency(agg$Annual), 1)
})

test_that("tsaggregates with only the original period returns the series", {
  agg <- tsaggregates(USAccDeaths, aggregatelist = 1)
  expect_length(agg, 1)
  expect_equal(agg[[1]], USAccDeaths)
})
