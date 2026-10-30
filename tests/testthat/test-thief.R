test_that("thief returns a forecast object of the right shape", {
  z <- thief(USAccDeaths, h = 12)
  expect_s3_class(z, "forecast")
  expect_equal(z$method, "THieF-ETS")
  expect_length(z$mean, 12)
  expect_equal(tsp(z$mean)[1], tsp(USAccDeaths)[2] + 1 / 12)
  expect_equal(frequency(z$mean), 12)
  expect_equal(z$x, USAccDeaths)
  expect_length(z$fitted, length(USAccDeaths))
  expect_equal(tsp(z$fitted), tsp(USAccDeaths))
  expect_false(anyNA(z$mean))
})

test_that("thief default horizon is two seasonal periods", {
  expect_length(thief(USAccDeaths, usemodel = "naive")$mean, 24)
})

test_that("thief works with each model and combination", {
  for (model in c("arima", "theta", "naive", "snaive")) {
    z <- thief(USAccDeaths, h = 6, usemodel = model)
    expect_equal(z$method, paste0("THieF-", toupper(model)))
    expect_length(z$mean, 6)
    expect_false(anyNA(z$mean))
  }
  for (comb in c("struc", "mse", "ols", "bu")) {
    z <- thief(USAccDeaths, h = 6, comb = comb, usemodel = "naive")
    expect_length(z$mean, 6)
    expect_false(anyNA(z$mean))
  }
})

test_that("thief works with GLS combinations given enough data", {
  set.seed(1)
  y <- ts(rnorm(240, 100), frequency = 12)
  # naive and snaive residuals contain NAs at the start
  for (model in c("ets", "naive", "snaive")) {
    for (comb in c("shr", "sam")) {
      z <- suppressWarnings(thief(y, h = 6, comb = comb, usemodel = model))
      expect_length(z$mean, 6)
      expect_false(anyNA(z$mean))
    }
  }
})

test_that("thief with bottom-up combination on naive forecasts equals naive", {
  z <- thief(USAccDeaths, h = 12, comb = "bu", usemodel = "naive")
  expect_equal(as.numeric(z$mean), rep(as.numeric(USAccDeaths[72]), 12))
})

test_that("thief accepts a user-supplied forecast function and aggregatelist", {
  fnaive <- function(y, h, ...) forecast::naive(y, h = h, ...)
  z <- thief(USAccDeaths, h = 12, forecastfunction = fnaive)
  expect_s3_class(z, "forecast")
  expect_equal(z$method, "THieF-FNAIVE")
  z2 <- thief(USAccDeaths, h = 12, usemodel = "naive", aggregatelist = c(1, 3, 12))
  expect_length(z2$mean, 12)
})

test_that("thief handles AEdemand data", {
  z <- thief(AEdemand[, 1], h = 12, usemodel = "naive")
  expect_length(z$mean, 12)
})

test_that("thief validates input", {
  expect_error(thief(1:48), "must be a time series")
  expect_error(thief(AEdemand[, 1:2]), "univariate")
  expect_error(thief(ts(1:48)), "greater than 1")
  expect_error(thief(ts(1:10, frequency = 12)), "at least 2 periods")
  expect_error(thief(USAccDeaths, aggregatelist = c(2, 4)), "aggregatelist must include")
  expect_error(thief(USAccDeaths, comb = "foo"), "should be one of")
  expect_error(thief(USAccDeaths, usemodel = "foo"), "should be one of")
})

test_that("the m argument of thief is deprecated", {
  expect_warning(
    z <- thief(USAccDeaths, m = 12, h = 6, usemodel = "naive"),
    "deprecated"
  )
  expect_equal(z$mean, thief(USAccDeaths, h = 6, usemodel = "naive")$mean)
  expect_error(suppressWarnings(thief(USAccDeaths, m = 4)), "m must equal frequency")
})

test_that("thief passes ... to theta, naive and snaive", {
  for (model in c("theta", "naive", "snaive")) {
    plain <- thief(USAccDeaths, h = 6, usemodel = model)
    boxcox <- thief(USAccDeaths, h = 6, usemodel = model, lambda = 0, biasadj = TRUE)
    expect_false(isTRUE(all.equal(plain$mean, boxcox$mean)))
  }
})

test_that("thief handles a forecastfunction returning a single-level lower vector", {
  f80 <- function(y, h, ...) {
    fc <- forecast::naive(y, h = h, level = 80)
    fc$lower <- as.numeric(fc$lower)
    fc
  }
  expect_no_error(z <- thief(USAccDeaths, h = 12, forecastfunction = f80))
  expect_length(z$mean, 12)
})
