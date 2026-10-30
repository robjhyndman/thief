test_that("the AEdemand example runs", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  pdf(NULL)
  on.exit(dev.off())

  p <- ggplot2::autoplot(AEdemand, xlab = "Year", ylab = "Demand ('000)")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  total <- AEdemand[, 12]
  totalagg <- tsaggregates(total)
  expect_named(totalagg, c("Weekly", "2-Weekly", "4-Weekly", "Quarterly", "Biannual", "Annual"))
  expect_no_error(plot(totalagg, main = "Total Emergency Admissions"))

  base <- list()
  for (i in 1:5) {
    base[[i]] <- forecast::forecast(forecast::auto.arima(totalagg[[i]]))
  }
  base[[6]] <- forecast::forecast(forecast::auto.arima(totalagg[[6]]), h = 2)
  reconciled <- reconcilethief(base)
  expect_equal(unname(vapply(reconciled, function(f) frequency(f$mean), numeric(1))), c(52, 26, 13, 4, 2, 1))
  expect_equal(sum(reconciled[[1]]$mean), sum(reconciled[[6]]$mean))
})
