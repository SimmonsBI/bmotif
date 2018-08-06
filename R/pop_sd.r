pop_sd <- function(v) {
  # compute population sd instead of sample sd
  # built in sd uses denominator n - 1
  n <- length(v)
  return(stats::sd(v) * sqrt((n - 1) / n))
}
