test_that("function can assess if distributions are equal", {
  # Create a dataframe that should pass
  pass <- data.frame(
    group = c(rep("A", 50), rep("B", 50)),
    value = c(rep(1, 100))
  )

  fail <- data.frame(
    group = c(rep("A", 50), rep("B", 50)),
    value = c(rep(1, 50), rep(10, 50))
  )

  expect_true(
    verify_equal_disitributions(
      pass,
      grouping_colname = "group",
      metric_colname = "value"
    )
  )

  expect_false(
    verify_equal_disitributions(
      fail,
      grouping_colname = "group",
      metric_colname = "value"
    )
  )
})
