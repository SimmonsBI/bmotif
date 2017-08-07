mprofile <- function(M, normalise){
  #' Count bipartite motifs
  #'
  #' Counts occurrences of motifs in a bipartite network
  #' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
  #' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes do interact, and 0
  #' otherwise. Formally, M is an incidence matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
  #' If interactions are weighted (matrix elements are greater than 1), the function will automatically convert the matrix to a binary
  #' matrix.
  #' @param normalise Logical; should motif frequencies be normalised to control for network size?
  #' @details Counts the number of times each of the 17 bipartite motifs up to five nodes occurs in a network.
  #' @return
  #' By default (\code{normalise} = FALSE), \code{mprofile} returns a data frame with 17 rows (one for each motif) and 3 columns.
  #' The first column (\code{motif}) indicates the motif ID as described in Simmons et al. (2017) (and originally in Appendix 1 of Baker et al. (2015)). The second column
  #' (\code{nodes}) indicates how many nodes the motif contains. The third column (\code{frequency}) is the number of times each motif appears in the network.
  #'
  #' If \code{normalise} = TRUE, three additional columns are added to the output data frame, each corresponding to a different method of normalising motif
  #' frequencies. The first column (\code{normalise_sum}) converts each frequency to a relative frequency by expressing counts as a proportion of the total number
  #' of motifs in the network. The second column (\code{normalise_sizeclass}) uses a similar approach, but expresses counts as a proportion of the total number of
  #' motifs within each motif size class (the number of nodes a motif contains). For example, the relative frequency of all two-node motifs will sum to one,
  #' as will the relative frequency of all three-, four- and five- node motifs. The final column (\code{normalise_nodesets}) expresses frequencies as the number
  #' of node sets that are involved in a motif as a proportion of the number of node sets that could be involved in that motif. For example, in a motif
  #' with three nodes in one level (A) and two nodes in the other level (P), the maximum number of node sets which could be involved in the motif is
  #' given by the product of binomial coefficients, choosing three nodes from A and two from P.
  #' @export
  #' @references
  #' Baker, N., Kaartinen, R., Roslin, T., and Stouffer, D. B. (2015). Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2):130–139.
  #'
  #' Simmons, B I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and di Clemente, R. bmotif: a package for counting motifs in bipartite networks
  #' @examples
  #' set.seed(123)
  #' row <- 100
  #' col <- 100
  #' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
  #' mprofile(M = m, normalise = TRUE)


  # check inputs
  if(class(M) != "matrix"){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length zero")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(class(normalise) != "logical"){stop("'normalise' must be of class 'logical' i.e. TRUE or FALSE")} # make sure normalise is logical i.e. TRUE or FALSE

  # clean matrix
  M[M > 0] <- 1 # ensure M is binary
  dimnames(M) <- NULL # strip row and column names

  # calculate inputs
  z <- dim(M)[1]
  p <- dim(M)[2]
  Tz <- M %*% t(M)
  Tp <- t(M) %*% M
  lTz <- dim(Tz)[2]
  lTp <- dim(Tp)[2]

  # create results container
  if(normalise == TRUE){
    out <- data.frame(motif = 1:17, nodes = c(1,rep(2,2),rep(4,4),rep(5,10)), frequency = NA, normalise_sum = NA, normalise_sizeclass = NA, normalise_nodesets = NA)
  } else {
    out <- data.frame(motif = 1:17, nodes = c(1,rep(2,2),rep(4,4),rep(5,10)), frequency = NA)
  }

  # count motifs
  for(i in 1:17){
    out[i,"frequency"] <- countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp)
  }

  # normalisations
  if(normalise == TRUE){
    # calculate normalised frequency across all motifs
    out$normalise_sum <- out$frequency/sum(out$frequency)

    # calculate normalised frequency within a motif size class
    out$normalise_sizeclass <- do.call("c", lapply(split(out, out$nodes),
                                                   function(df) df$frequency/sum(df$frequency))
    )

    # calculate normalised frequency as proportion of possible node sets
    sets <- node_sets(M)
    out$normalise_nodesets <- out$frequency/sets
  }

  # output
  return(out)
}
