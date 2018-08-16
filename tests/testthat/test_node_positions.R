context("node_positions")

test_that("In a motif, each node is only in one position",{
  mp <- lapply(1:44, function(x) motif_info(x, link = FALSE)) # motif positions
  for(i in 1:length(mp)){
    ps <- mp[[i]]
    mot <- motifs[[i]]
    nps <- node_positions(M = mot, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    nps <- nps[,ps]
    expect_identical(all(rowSums(nps) == 1), TRUE)
  }

  # tests with block matrices
  mlist2 <- motifs[1]
  bm2 <- block_matrix(mlist2)
  np <- node_positions(M = bm2, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,1:2]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist3 <- motifs[2:3]
  bm3 <- block_matrix(mlist3)
  np <- node_positions(M = bm3, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,3:6]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist4 <- motifs[4:7]
  bm4 <- block_matrix(mlist4)
  np <- node_positions(M = bm4, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,7:16]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist5 <- motifs[8:17]
  bm5 <- block_matrix(mlist5)
  np <- node_positions(M = bm5, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,17:46]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist6 <- motifs[18:44]
  bm6 <- block_matrix(mlist6)
  np <- node_positions(M = bm6, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,47:148]
  expect_identical(all(rowSums(np) == 1), TRUE)
})

test_that("Error messages work",{
  expect_error(object = node_positions("a"), "'M' must be an object of class 'matrix'")
  expect_error(object = node_positions(matrix("a",3,3)), "Elements of 'M' must be numeric")
  expect_error(object = node_positions(matrix(-1,3,3)), "Elements of 'M' must be greater than or equal to zero")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = 7, weights_method = "none", weights_combine = "none", normalisation = "none"), "'level' must be of class 'character'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "wrong", weights_method = "none", weights_combine = "none", normalisation = "none"), "'level' must equal 'rows', 'columns' or 'all'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = 7), "'normalisation' must be of class 'character'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "wrong"), "'normalisation' must equal 'none','sum','size class', 'size class_plus1', 'size class_NAzero', 'position','levelsize','levelsize_plus1', 'levelsize_NAzero', 'motif','motif_plus1','motif_NAzero'")
})

test_that("Matches validated results",{
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none"),positions_T_all_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_all_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_T_all_sc)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_T_rows_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_rows_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_T_rows_sc)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_T_columns_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_columns_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_T_columns_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_all_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_all_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_F_all_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_rows_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_rows_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_F_rows_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_columns_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_columns_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "size class"), positions_F_columns_sc)
})
