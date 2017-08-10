multtensor <- function(A, B){ # Cikl = Aij * Bjkl
  #' @importFrom tensor tensor
  tensor(A,B,2,1)
}
