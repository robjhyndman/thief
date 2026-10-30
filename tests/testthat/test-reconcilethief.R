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

make_fc <- function(y = USAccDeaths, h = 12) {
  lapply(tsaggregates(y), function(x) {
    forecast::forecast(forecast::ets(x), h = h / (12 / frequency(x)))
  })
}

test_that("reconcilethief makes forecasts coherent across all combinations", {
  fc <- make_fc()
  for (comb in c("struc", "mse", "ols")) {
    mse <- vapply(fc, function(f) mean(residuals(f)^2), numeric(1))
    rec <- reconcilethief(fc, comb = comb, mse = mse)
    expect_length(rec, length(fc))
    expect_s3_class(rec[[1]], "forecast")
    # Aggregates equal the sums of the reconciled bottom level
    expect_equal(as.numeric(rec[[length(rec)]]$mean), sum(rec[[1]]$mean))
    expect_equal(sum(rec[[3]]$mean), sum(rec[[1]]$mean))
  }
})

test_that("reconcilethief GLS combinations give coherent forecasts", {
  set.seed(1)
  fc <- make_fc(ts(rnorm(240, 100), frequency = 12))
  for (comb in c("shr", "sam")) {
    rec <- suppressWarnings(reconcilethief(fc, comb = comb))
    expect_equal(sum(rec[[3]]$mean), sum(rec[[1]]$mean))
  }
})

test_that("reconcilethief bottom-up keeps the bottom forecasts", {
  fc <- make_fc()
  rec <- reconcilethief(fc, comb = "bu")
  expect_equal(as.numeric(rec[[1]]$mean), as.numeric(fc[[1]]$mean))
})

test_that("reconcilethief returnall = FALSE returns only the bottom level", {
  fc <- make_fc()
  rec <- reconcilethief(fc, returnall = FALSE)
  expect_s3_class(rec, "forecast")
  expect_length(rec$mean, 12)
})

test_that("reconcilethief accepts plain ts lists", {
  fc <- lapply(make_fc(), function(f) f$mean)
  rec <- reconcilethief(fc, comb = "ols")
  expect_s3_class(rec, "tsaggregates")
  expect_equal(sum(rec[[length(rec)]]), sum(rec[[1]]))
  bts <- reconcilethief(fc, comb = "ols", returnall = FALSE)
  expect_s3_class(bts, "ts")
  expect_length(bts, 12)
})

test_that("reconcilethief validates input", {
  fc <- lapply(make_fc(), function(f) f$mean)
  expect_error(reconcilethief(fc[-length(fc)]), "Minimum seasonal period")
  bad <- fc
  bad[[1]] <- bad[[1]][1:6]
  expect_error(reconcilethief(bad), "same equivalent lengths")
  expect_error(reconcilethief(fc, comb = "shr"), "GLS needs residuals")
  expect_error(reconcilethief(fc, comb = "foo"), "should be one of")
})

test_that("reconcilethief does not depend on the order of the input list", {
  fc <- make_fc()
  mse <- vapply(fc, function(f) mean(residuals(f)^2), numeric(1))
  for (comb in c("struc", "mse")) {
    expect_no_warning(rec <- reconcilethief(fc, comb = comb, mse = mse))
    expect_no_warning(rec2 <- reconcilethief(rev(fc), comb = comb, mse = rev(mse)))
    expect_equal(rec2, rec)
  }
})

test_that("reconcilethief gives a clear error when comb = 'mse' has no mse", {
  expect_error(reconcilethief(make_fc(), comb = "mse"), "needs the mse argument")
})
