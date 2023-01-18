test_that("verify_normality correctly assesses distributions", {
  # For Petal.Width in iris dataset, virginica produces a non-significant
  # result, however setosa is incredibly significant.
  # For Sepal.With, both virginica and versicolor produce non-significant
  # results.

  # Create dataframe that will pass
  pass <- iris %>%
    dplyr::filter(Species %in% c("versicolor", "virginica"))

  # Create dataframe that will fail
  fail <- iris %>%
    dplyr::filter(Species %in% c("versicolor", "setosa"))

  expect_true(
    verify_normality(
      pass,
      grouping_colname = "Species",
      metric_colname = "Sepal.Width",
      subset = "versicolor"
    )
  )

  expect_false(
    verify_normality(
      fail,
      grouping_colname = "Species",
      metric_colname = "Petal.Width",
      subset = "versicolor"
    )
  )
})
