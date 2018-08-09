context("mcount")

test_that("mcount behaves as expected",{
  expect_error(object = mcount("a"), "'M' must be an object of class 'matrix'")
  expect_error(object = mcount(matrix("a",3,3)), "Elements of 'M' must be numeric")
  expect_error(object = mcount(matrix(-1,3,3)), "Elements of 'M' must be greater than or equal to zero")
  expect_error(object = mcount(matrix(1,3,3), six_node = TRUE, normalisation = 7), "'normalisation' must be of class 'logical' i.e. TRUE or FALSE")
  expect_error(object = mcount(matrix(1,3,3), six_node = 7, normalisation = TRUE), "'six_node' must be of class 'logical' i.e. TRUE or FALSE")

  expect_equal(object = mcount(M = mat, six_node = FALSE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_SF_NF)
  expect_equal(object = mcount(M = mat, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_ST_NF)
  expect_equal(object = mcount(M = mat, six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE), mcount_SF_NT)
  expect_equal(object = mcount(M = mat, six_node = TRUE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE), mcount_ST_NT)
  expect_equal(mcount(M = t(mat), six_node = TRUE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_ST_NT)
  expect_equal(mcount(M = t(mat), six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_SF_NT)
  expect_equal(mcount(M = t(mat), six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_ST_NF)
  expect_equal(mcount(M = t(mat), six_node = FALSE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_SF_NF)
})
