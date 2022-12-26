

test_that("group column is assigned correctly to subset and whole without subset", {
  test_df <- data.frame(
    subtype = c("A", "A", "B", "B", "C", "C"),
    metric = c(1, 1, 2, 2, 3, 3)
  )

  expected_result <- test_df %>%
    dplyr::mutate(group = c(rep("subset", 2), rep("whole_without_subset", 4)))

  expect_equal(
    define_subset_against_whole(test_df, "subtype", "A", "metric"),
    expected_result
  )
})
