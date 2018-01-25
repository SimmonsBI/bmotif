test_that("mcount behaves as expected",{
  expect_identical(object = mcount(M = mat, six_node = FALSE, normalise = FALSE), mcount_SF_NF)
  expect_identical(object = mcount(M = mat, six_node = TRUE, normalise = FALSE), mcount_ST_NF)
  expect_identical(object = mcount(M = mat, six_node = FALSE, normalise = TRUE), mcount_SF_NT)
  expect_identical(object = mcount(M = mat, six_node = TRUE, normalise = TRUE), mcount_ST_NT)
})
