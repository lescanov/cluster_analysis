test_that("function can correctly assess equality of variances", {
  # Turns out that in the iris dataset, Species yields
  # equal variances and Sepal.Lenggth is unequal variances
  # as per levene_test from rstatix package
  expect_true(
    verify_equality_of_variances(iris, "Species", "Sepal.Width")
  )

  expect_false(
    verify_equality_of_variances(iris, "Species", "Sepal.Length")
  )
})
