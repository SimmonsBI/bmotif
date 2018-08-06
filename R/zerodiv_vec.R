zerodiv_vec <- function(v,w) {
  a <- cbind(v,w)
  apply(a, 1, zerodiv)
}