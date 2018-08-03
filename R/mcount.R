mcount <- function(M, six_node, normalisation){
  #' Count bipartite motifs
  #'
  #' Counts occurrences of motifs in a bipartite network
  #' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
  #' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes do interact, and 0
  #' otherwise. Formally, M is an incidence matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
  #' If interactions are weighted (non-zero matrix elements take values other than 1), the function will automatically convert the matrix to a binary
  #' matrix.
  #' @param six_node Logical; should six node motifs be counted?
  #' @param normalisation Logical; should motif frequencies be normalised to control for network size?
  #' @details Counts the number of times each of the 17 motifs up to five nodes (if \code{six_node} = FALSE), or 44 motifs up to six nodes (if \code{six_node} = TRUE), occurs in a network.
  #'
  #' Larger networks tend to contain more motifs. Controlling for this effect by normalising motif counts is important if different sized networks are being compared.
  #' If \code{normalisation} = TRUE, motif frequencies are normalised in three ways. The first method ("normalise_sum") converts each frequency to a relative frequency by expressing counts as
  #' a proportion of the total number of motifs in the network. The second method ("normalise_sizeclass") uses a similar approach, but expresses counts as a proportion of the total number of
  #' motifs within each motif size class (the number of nodes a motif contains). For example, the relative frequency of all two-node motifs will sum to one,
  #' as will the relative frequency of all three-, four-, five- and six-node motifs. The final method ("normalise_nodesets") expresses frequencies as the number
  #' of node sets that are involved in a motif as a proportion of the number of node sets that could be involved in that motif (Poisot and Stouffer, 2017). For example, in a motif
  #' with three nodes in one level (A) and two nodes in the other level (P), the maximum number of node sets which could be involved in the motif is
  #' given by the product of binomial coefficients, choosing three nodes from A and two from P.
  #'
  #' Warning: including six node motifs is fine for most networks. However, for large networks, counting six node motifs can be slow and memory intensive. In some cases, R can crash if there is not enough memory.
  #' @return
  #' Returns a data frame with one row for each motif: either 17 rows (if \code{six_node} = FALSE) or 44 rows (if \code{six_node} = TRUE). The data frame has three columns.
  #' The first column ("motif") indicates the motif ID as described in Simmons et al. (2017) (and originally in Appendix 1 of Baker et al. (2015)). The second column
  #' ("nodes") indicates how many nodes the motif contains. The third column ("frequency") is the number of times each motif appears in the network.
  #'
  #' If \code{normalisation} = TRUE, three additional columns are added to the output data frame, each corresponding to a different method of normalising motif
  #' frequencies as described above.
  #' @export
  #' @references
  #' Baker, N., Kaartinen, R., Roslin, T., and Stouffer, D. B. (2015). Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2):130–139.
  #'
  #' Poisot, T. & Stouffer, D. (2016). How ecological networks evolve. bioRxiv.
  #'
  #' Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356
  #' @examples
  #' set.seed(123)
  #' row <- 10
  #' col <- 10
  #' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
  #' mcount(M = m, six_node = TRUE, normalisation = TRUE)

  # check inputs
  if(class(M) != "matrix"){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length zero")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(class(normalisation) != "logical"){stop("'normalisation' must be of class 'logical' i.e. TRUE or FALSE")} # make sure normalisation is logical i.e. TRUE or FALSE
  if(class(six_node) != "logical"){stop("'six_node' must be of class 'logical' i.e. TRUE or FALSE")} # make sure six_node is logical i.e. TRUE or FALSE

  # clean matrix
  M[M > 0] <- 1 # ensure M is binary
  dimnames(M) <- NULL # strip row and column names

  # calculate inputs
  p <- dim(M)[2]
  z <- dim(M)[1]
  J <- matrix(rep(1, z * p), nrow = z, ncol = p)
  JP <- matrix(rep(1, p * p), nrow = p, ncol = p)
  JZ <- matrix(rep(1, z * z), nrow = z, ncol = z)
  MT <- t(M)
  N <- J - M
  NT <- t(N)
  P <- MT %*% M
  Q <- MT %*% N
  R <- NT %*% M
  Z <- M %*% MT
  Y <- M %*% NT
  X <- N %*% MT
  dP <- apply(M, MARGIN = 2, sum)
  jP <- rep(1, p)
  dZ <- apply(M, MARGIN = 1, sum)
  jZ <- rep(1, z)

  if(six_node == TRUE){
    if (p < z) {
      J3 <- array(rep(1, p * p * p), c(p, p, p))
      AP <- maketensor(M, M)
      BP <- maketensor(M, N)
      CP <- maketensor(N, M)
      DP <- maketensor(N, N)
      MA <- tensor::tensor(MT, AP, 2, 1)
      MB <- tensor::tensor(MT, BP, 2, 1)
      MC <- tensor::tensor(MT, CP, 2, 1)
      MD <- tensor::tensor(MT, DP, 2, 1)
      Na <- tensor::tensor(NT, AP, 2, 1)
      NB <- tensor::tensor(NT, BP, 2, 1)
      NC <- tensor::tensor(NT, CP, 2, 1)
      K3 <- J3
      for (i in 1 : p){
        for (j in 1 : p){
          K3[i,j,j] <- 0
          K3[j,i,j] <- 0
          K3[j,j,i] <- 0
        }
      }
    }
    if (p >= z) {
      J3 <- array(rep(1, z * z * z), c(z, z, z))
      AP <- maketensor(MT, MT)
      BP <- maketensor(MT, NT)
      CP <- maketensor(NT, MT)
      DP <- maketensor(NT, NT)
      MA <- tensor::tensor(M, AP, 2, 1)
      MB <- tensor::tensor(M, BP, 2, 1)
      MC <- tensor::tensor(M, CP, 2, 1)
      MD <- tensor::tensor(M, DP, 2, 1)
      Na <- tensor::tensor(N, AP, 2, 1)
      NB <- tensor::tensor(N, BP, 2, 1)
      NC <- tensor::tensor(N, CP, 2, 1)
      K3 <- J3
      for (i in 1 : z){
        for (j in 1 : z){
          K3[i,j,j] <- 0
          K3[j,i,j] <- 0
          K3[j,j,i] <- 0
        }
      }
    }
    MA <- MA * K3
    MB <- MB * K3
    MC <- MC * K3
    MD <- MD * K3
    Na <- Na * K3
    NB <- NB * K3
    NC <- NC * K3
  }

  # create results container
  if(six_node == FALSE){
    if(normalisation == TRUE){
      out <- data.frame(motif = 1:17, nodes = c(2,rep(3,2),rep(4,4),rep(5,10)), frequency = NA, normalise_sum = NA, normalise_sizeclass = NA, normalise_nodesets = NA)
    } else {
      out <- data.frame(motif = 1:17, nodes = c(2,rep(3,2),rep(4,4),rep(5,10)), frequency = NA)
    }
  } else {
    if(normalisation == TRUE){
      out <- data.frame(motif = 1:44, nodes = c(2,rep(3,2),rep(4,4),rep(5,10),rep(6,27)), frequency = NA, normalise_sum = NA, normalise_sizeclass = NA, normalise_nodesets = NA)
    } else {
      out <- data.frame(motif = 1:44, nodes = c(2,rep(3,2),rep(4,4),rep(5,10),rep(6,27)), frequency = NA)
    }
  }

  # count motifs
  if(six_node == FALSE){
    for(i in 1:17){
      out[i,"frequency"] <- countmotif(x = M, motif =  i, z = z, p = p, JP = JP, JZ = JZ, P = P, Q = Q, R = R, Z = Z, Y = Y, X = X, dP = dP, jP = jP, dZ = dZ, jZ = jZ)
    }
  } else {
    for(i in 1:44){
      out[i,"frequency"] <- countmotif(x = M, motif =  i, z = z, p = p, JP = JP, JZ = JZ, P = P, Q = Q, R = R, Z = Z, Y = Y, X = X, dP = dP, jP = jP, dZ = dZ, jZ = jZ, J3 = J3, MA = MA, MB = MB, MC = MC, MD = MD, Na = Na, NB = NB, NC = NC)
    }
  }

  # normalisations
  if(normalisation == TRUE){
    # calculate normalised frequency across all motifs
    out$normalise_sum <- out$frequency/sum(out$frequency)

    # calculate normalised frequency within a motif size class
    out$normalise_sizeclass <- do.call("c", lapply(split(out, out$nodes),
                                                   function(df) df$frequency/sum(df$frequency))
    )

    # calculate normalised frequency as proportion of possible node sets
    sets <- node_sets(M, six_node = six_node)
    out$normalise_nodesets <- out$frequency/sets
  }

  # output
  return(out)
}
