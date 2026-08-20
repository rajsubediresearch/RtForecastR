test_that("rt_forecast runs end-to-end on the bundled example data", {
  data(measles_cdmx)
  fit <- rt_forecast(measles_cdmx$time, measles_cdmx$cases,
                      mean_GI = 11/7, var_GI = (4/7)^2)
  expect_s3_class(fit, "rtforecast")
  expect_true(nrow(fit$results) > 0)
  expect_true(nrow(fit$predictions) > 0)
  expect_equal(nrow(fit$forecast), 1)
  expect_true(fit$forecast$time > max(measles_cdmx$time))
  expect_true(fit$forecast$lo95 <= fit$forecast$pred_mean)
  expect_true(fit$forecast$hi95 >= fit$forecast$pred_mean)
})

test_that("rt_forecast handles a larger-incidence series via maxI", {
  data(jalisco)
  fit <- rt_forecast(jalisco$time, jalisco$cases, mean_GI = 11/7, var_GI = (4/7)^2)
  expect_s3_class(fit, "rtforecast")
  expect_true(!is.null(fit$forecast_quantiles))
})

test_that("forecast_quantiles exactly matches fit$forecast at the 95% and 50% levels", {
  # Regression test for a bug where forecast_quantiles was computed from
  # a disconnected in-sample calculation rather than the same
  # distribution used for fit$forecast itself.
  data(measles_cdmx)
  fit <- rt_forecast(measles_cdmx$time, measles_cdmx$cases,
                      mean_GI = 11/7, var_GI = (4/7)^2)
  expect_equal(as.numeric(fit$forecast_quantiles[["0.025"]]),
               c(fit$forecast$lo95, fit$forecast$hi95))
  expect_equal(as.numeric(fit$forecast_quantiles[["0.25"]]),
               c(fit$forecast$lo50, fit$forecast$hi50))
  # no duplicate levels
  expect_equal(length(names(fit$forecast_quantiles)),
               length(unique(names(fit$forecast_quantiles))))
})

test_that("rt_forecast errors on mismatched input lengths", {
  expect_error(rt_forecast(1:5, 1:4, mean_GI = 1, var_GI = 1))
})

test_that("compute_lambda's first element is NA and length matches input", {
  lam <- compute_lambda(c(1, 5, 10, 8), mean_GI = 1.5, var_GI = 0.5)
  expect_true(is.na(lam[1]))
  expect_equal(length(lam), 4)
})
