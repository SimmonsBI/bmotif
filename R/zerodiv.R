zerodiv <- function(vec) {
  x <- vec[1]
  y <- vec[2]
  ifelse(y == 0, return (0), return(x / y))
}