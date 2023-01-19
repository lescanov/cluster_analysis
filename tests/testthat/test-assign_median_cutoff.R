test_that("each sample is correctly assessed as above or below median", {
  # Creating test dataframe
  test <- data.frame(
    sample = c(1, 2, 3, 4, 5),
    value = c(1, 2, 3, 4, 5)
  )

  expected_result <- test %>%
    mutate(
      median_status = c("below", "below", "below", "above", "above")
    )

  expect_equal(
    assign_median_cutoff(test, "value"),
    expected_result
  )
})
