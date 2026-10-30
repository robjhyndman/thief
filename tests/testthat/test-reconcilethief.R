test_that("reconcilethief handles length-1 series with rounding differences in start (#19)", {
  set.seed(1)
  # Start times of tsaggregates() and forecast() differ by ~1e-13 here
  y <- ts(rnorm(101) + 10, start = c(1990, 1), frequency = 12)
  fc <- lapply(tsaggregates(y), function(x) {
    forecast::forecast(forecast::ets(x), h = frequency(x))
  })
  expect_false(identical(tsp(fc[[6]]$mean)[1], tsp(tsaggregates(fc[[1]]$mean)[[6]])[1]))
  expect_no_warning(rec <- reconcilethief(fc))
  annual <- rec[[length(rec)]]
  expect_length(annual$mean, 1)
  expect_false(anyNA(annual$lower))
  expect_false(anyNA(annual$upper))
  expect_equal(dim(annual$lower), c(1, 2))
  # Intervals are shifted by the same adjustment as the mean
  adj <- as.numeric(annual$mean) - as.numeric(fc[[6]]$mean)
  expect_equal(as.numeric(annual$lower), as.numeric(fc[[6]]$lower) + adj)
})
