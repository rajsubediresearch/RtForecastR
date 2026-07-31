test_that("mae/mse/rmse compute correctly", {
  expect_equal(mae(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_equal(mae(c(0, 0), c(3, 4)), 3.5)
  expect_equal(mse(c(0, 0), c(3, 4)), (9 + 16) / 2)
  expect_equal(rmse(c(0, 0), c(3, 4)), sqrt((9 + 16) / 2))
})

test_that("mae/mse error on length mismatch", {
  expect_error(mae(1:3, 1:2))
  expect_error(mse(1:3, 1:2))
})

test_that("mae ignores NA pairs", {
  expect_equal(mae(c(1, 2, NA), c(1, 2, 5)), 0)
})

test_that("coverage computes the correct fraction", {
  expect_equal(coverage(c(5, 5, 5), c(0, 0, 0), c(10, 10, 10)), 1)
  expect_equal(coverage(c(50, 50, 50), c(0, 0, 0), c(10, 10, 10)), 0)
  expect_equal(coverage(c(5, 50), c(0, 0), c(10, 10)), 0.5)
})

test_that("interval_score matches the closed-form for inside/outside points", {
  # observed inside interval: score == width only
  expect_equal(interval_score(5, 2, 8, 0.05), 8 - 2)
  # observed below interval: width + penalty
  expect_equal(interval_score(20, 2, 8, 0.05), (8 - 2) + (2/0.05) * (20 - 8))
})

test_that("wis returns per-point scores and a mean_wis attribute", {
  q <- list("0.025" = matrix(c(2, 8, 15, 25), nrow = 2),
            "0.25"  = matrix(c(4, 6, 18, 22), nrow = 2))
  result <- wis(observed = c(5, 20), median_pred = c(5, 21), quantiles = q)
  expect_equal(nrow(result), 2)
  expect_true(all(result$wis >= 0))
  expect_true(is.numeric(attr(result, "mean_wis")))
})
