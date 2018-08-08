context("Weighted link positions")

# Define some functions
rbm <- function (m, n, one_prob = 0.5) {
  matrix(sample(0:1, m * n, replace = TRUE, prob = c(1- one_prob,one_prob)), m, n)
}

# Tests
test_that("Check binary matrices", {
  for (i in 1:10) {
    W <- rbm(20,20)
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    expect(all(lw == lp))
  }
})

test_that("Check matrices with equal weights", {
  for (i in 1:5) {
    W <- 3 * rbm(20,20)
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    expect(all(lw == 3 * lp))
  }
  for (i in 1:5) {
    W <- 0.5 * rbm(20,20)
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    expect(all(lw == 0.5 * lp))
  }
})


test_that("Check matrices with one weighted link", {
  for (i in 1:5) {
    W <- rbm(20,20)
    W[1,1] <- 10
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    expect(all(lw[1,] == 10 * lp[1,]))
    expect(all(lw[-1,] == lp[-1,]))
  }
  for (i in 1:5) {
    W <- rbm(20,20)
    W[3,3] <- 0.5
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    ind <- which(rownames(lp) == "r3 -- c3")
    expect(all(lw[ind,] == 0.5 * lp[ind,]))
    expect(all(lw[-ind,] == lp[-ind,]))
  }
})

test_that("Check matrices with several weighted links", {
  for (i in 1:5) {
    W <- rbm(20,20)
    W[7,8] <- 1.5
    W[4,6] <- 0.7
    lp <- link_positions(W, weights = FALSE)
    lw <- link_positions(W, weights = TRUE)
    ind <- which(rownames(lp) == "r7 -- c8")
    ind2 <- which(rownames(lp) == "r4 -- c6")
    expect(all(lw[ind,] == 1.5 * lp[ind,]))
    expect(all(lw[ind2,] == 0.7 * lp[ind2,]))
    expect(all(lw[-c(ind, ind2),] == lp[-c(ind, ind2),]))
  }
})
