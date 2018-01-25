test_that("mcount behaves as expected",{
  expect_identical(object = mcount(M = mat, six_node = FALSE, normalise = FALSE), mcount_SF_NF)
  expect_identical(object = mcount(M = mat, six_node = TRUE, normalise = FALSE), mcount_ST_NF)
  expect_identical(object = mcount(M = mat, six_node = FALSE, normalise = TRUE), mcount_SF_NT)
  expect_identical(object = mcount(M = mat, six_node = TRUE, normalise = TRUE), mcount_ST_NT)
  expect_identical(mcount(M = t(mat), six_node = TRUE, normalise = TRUE), mcount_transpose_ST_NT)
  expect_identical(mcount(M = t(mat), six_node = FALSE, normalise = TRUE), mcount_transpose_SF_NT)
  expect_identical(mcount(M = t(mat), six_node = TRUE, normalise = FALSE), mcount_transpose_ST_NF)
  expect_identical(mcount(M = t(mat), six_node = FALSE, normalise = FALSE), mcount_transpose_SF_NF)

  expect_error(object = mcount("a"), "'M' must be an object of class 'matrix'")



  # if(class(M) != "matrix"){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  # if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  # if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length zero")} # make sure no elements of M have 0 length e.g. numeric(0)
  # if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  # if(class(normalise) != "logical"){stop("'normalise' must be of class 'logical' i.e. TRUE or FALSE")} # make sure normalise is logical i.e. TRUE or FALSE
  # if(class(six_node) != "logical"){stop("'six_node' must be of class 'logical' i.e. TRUE or FALSE")} # make sure six_node is logical i.e. TRUE or FALSE
})
