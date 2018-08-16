context("normalise_node_positions")

pc <- matrix(1, 10, 148, dimnames = list(c(paste0("r",1:5),paste0("c",1:5)),paste0("np",1:148)))

test_that("'none' doesn't change anything",{
  nnp <- normalise_node_positions(pc = pc, type = "none", six_node = TRUE)
  expect_identical(pc, nnp)
})

test_that("'sum' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sum", six_node = TRUE)
  expect_identical(all(rowSums(nnp) == 1), TRUE)
})

test_that("'position' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "position", six_node = TRUE)
  expect_identical(all(colSums(nnp) == 1), TRUE)
})

test_that("'sizeclass' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sizeclass", six_node = TRUE)
  expect_identical(all(rowSums(nnp[,1:2]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,3:6]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,7:16]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,17:46]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,47:148]) == 1), TRUE)
})

test_that("'levelsize' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "levelsize", six_node = TRUE)
  expect_identical(all(rowSums(nnp[,1:2]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,3:4]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,5:6]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,7:8]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,9:14]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,15:16]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,17:18]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,19:31]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,32:44]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,45:46]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,47:48]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,49:70]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,71:124]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,125:146]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,147:148]) == 1), TRUE)
})

test_that("'motif' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "motif", six_node = TRUE)
  mps <- lapply(1:44, function(x) motif_info(x, link = FALSE)) # motif positions
  for(i in mps){
    expect_identical(all(rowSums(nnp[,i]) == 1), TRUE)
  }
})

### ------ Testing _NAzero
test_that("'sizeclass_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(17,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "sizeclass", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "sizeclass_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

test_that("'levelsize_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "levelsize", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "levelsize_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

test_that("'motif_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "motif", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "motif_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

### ------ Testing _plus1
test_that("'sizeclass_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(17,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "sizeclass", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "sizeclass_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})

test_that("'levelsize_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "levelsize", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "levelsize_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})

test_that("'motif_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "motif", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "motif_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})
