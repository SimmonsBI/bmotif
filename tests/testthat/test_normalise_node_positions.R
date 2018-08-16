context("normalise_node_positions")

pc <- matrix(1, 10, 148, dimnames = list(c(paste0("r",1:5),paste0("c",1:5)),paste0("np",1:148)))

test_that("'none' doesn't change anything",{
  nnp <- normalise_node_positions(pc = pc, type = "none")
  expect_identical(pc, nnp)
})

test_that("'sum' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sum")
  expect_identical(all(rowSums(nnp) == 1), TRUE)
})

test_that("'position' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "position")
  expect_identical(all(colSums(nnp) == 1), TRUE)
})

test_that("'sizeclass' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sizeclass")
  expect_identical(all(rowSums(nnp[,1:2]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,3:6]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,7:16]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,17:46]) == 1), TRUE)
  expect_identical(all(rowSums(nnp[,47:148]) == 1), TRUE)
})

test_that("'levelsize' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "levelsize")
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
  nnp <- normalise_node_positions(pc = pc, type = "motif")
  mps <- lapply(1:44, function(x) motif_info(x, link = FALSE)) # motif positions
  for(i in mps){
    expect_identical(all(rowSums(nnp[,i]) == 1), TRUE)
  }
})
