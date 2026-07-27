# Tests for survival_utils.R

make_surv_df <- function() {
  set.seed(11)
  n <- 60
  data.frame(
    time = rexp(n, rate = 0.05) + 1,
    status = rbinom(n, 1, 0.6),
    age = rnorm(n, 60, 10),
    expr = rnorm(n),
    stringsAsFactors = FALSE
  )
}

test_that("validate_surv_df rejects bad input and accepts good input", {
  df <- make_surv_df()
  expect_invisible(validate_surv_df(df))

  bad_time <- df; bad_time$time[1] <- -5
  expect_error(validate_surv_df(bad_time), "> 0")

  bad_status <- df; bad_status$status[1] <- 2
  expect_error(validate_surv_df(bad_status), "0 \\(censored\\) or 1")

  expect_error(validate_surv_df(df[, c("time", "age")]), "must contain")
})

test_that("run_univariate_cox returns HR table for valid variables", {
  skip_if_not_installed("survival")
  df <- make_surv_df()
  out <- run_univariate_cox(df, c("age", "expr"))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_true(all(c("variable", "hazard_ratio", "HR_lower_95", "HR_upper_95", "pvalue") %in% colnames(out)))
  expect_true(all(out$hazard_ratio > 0))
  expect_true(all(out$HR_lower_95 <= out$HR_upper_95))
})

test_that("run_univariate_cox skips missing/too-small variables with warning", {
  skip_if_not_installed("survival")
  df <- make_surv_df()
  expect_warning(out <- run_univariate_cox(df, c("age", "nope")), "Variable not found")
  expect_equal(nrow(out), 1)  # only age
})

test_that("run_multivariate_cox fits joint model and validates vars", {
  skip_if_not_installed("survival")
  df <- make_surv_df()
  out <- run_multivariate_cox(df, c("age", "expr"))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_error(run_multivariate_cox(df, c("age", "nope")), "not found")
})

test_that("stratify_by_quantile makes balanced groups", {
  x <- 1:100
  g2 <- stratify_by_quantile(x, n_groups = 2)
  expect_equal(levels(g2), c("Low", "High"))
  expect_equal(as.numeric(table(g2)), c(50, 50))

  g4 <- stratify_by_quantile(x, n_groups = 4)
  expect_equal(length(levels(g4)), 4)
  expect_equal(as.numeric(table(g4)), rep(25, 4))

  expect_error(stratify_by_quantile(c(1, 1, 1), n_groups = 2), "Not enough unique")
})
