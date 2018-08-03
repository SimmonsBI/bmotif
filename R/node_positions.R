node_positions <- function(M, six_node, level = "all", normalisation = "none"){
  #' Calculate node position vectors
  #'
  #' Counts the frequency with which nodes occur in different positions within motifs.
  #' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
  #' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes do interact, and 0
  #' otherwise. Formally, M is an incidence matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
  #' If interactions are weighted (non-zero matrix elements take values other than 1), the function will automatically convert the matrix to a binary
  #' matrix.
  #' @param six_node Logical; should six node motifs be counted?
  #' @param  level Which node level should positions be calculated for: "rows", "columns" or "all"?  Defaults to "all".
  #' @param normalisation Which normalisation should be used: "none", "sum" or "size class"?  Defaults to "none".
  #' @details Counts the number of times each node in a network occurs in each of the 46 (if \code{six_node} = FALSE) or 148 (if \code{six_node} = TRUE) unique positions within motifs (to quantify a node's structural role).
  #'
  #' The \code{level} argument controls which node group positions are calculated for: "rows" returns position counts for all nodes in rows, "columns"
  #' returns counts for all nodes in columns, and "all" return counts for all nodes in the network.
  #'
  #' Nodes with more interactions will tend to appear in more positions. Normalisation helps control for this effect.
  #' "none" performs no normalisation and will return the raw position counts.
  #' "sum" divides position counts for each node by the total number of times that node appears in any position.
  #' "size class" divides position counts for each node by the total number of times that node appears in any position within the same motif size class.
  #'
  #' If a matrix is provided without row or column names, default names will be assigned: the first row will be called called 'r1', the second row will be called 'r2' and so on. Similarly, the first column will be called 'c1', the second column will be called 'c2' and so on.
  #'
  #' Warning: including six node motifs is fine for most networks. However, for large networks, counting six node motifs can be slow and memory intensive. In some cases, R can crash if there is not enough memory.
  #' @return
  #' Returns a data frame with one column for each node position: 46 columns if \code{six_node} is FALSE, and 148 columns if \code{six_node} is TRUE.
  #' Columns names are given as "npx" where x is the ID of the position as described in Simmons et al. (2017) (and originally in Appendix 1 of Baker et al. (2015))
  #'
  #' For a network with A rows and P columns, by default (where \code{level} = "all") the data frame has A + P rows, one for each node. If \code{level} = "rows", the data frame will have A rows, one for each row node;
  #' if \code{level} = "columns", it will have P rows, one for each column node.
  #'
  #' By default, the elements of this data frame will be the raw position counts. If \code{normalisation} is set to "sum" or "size class", the elements will be
  #' normalised position counts as described above.
  #' @export
  #' @references
  #' Baker, N., Kaartinen, R., Roslin, T., and Stouffer, D. B. (2015). Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2):130–139.
  #'
  #' Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356
  #' @examples
  #' set.seed(123)
  #' row <- 15
  #' col <- 15
  #' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
  #' rownames(m) <- paste0("r", 1:nrow(m)) # give the matrix row names
  #' colnames(m) <- paste0("c", 1:ncol(m)) # give the matrix column names
  #' node_positions(M = m, six_node = TRUE, level = "all", normalisation = "none")

  # check inputs
  if(class(M) != "matrix"){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length 0")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(class(level) != "character"){stop("'level' must be of class 'character'")} # make sure level is a character
  if(!level %in% c("rows","columns","all")){stop("'level' must equal 'rows', 'columns' or 'all'")} # make sure level equals 'rows', 'columns' or 'all'
  if(class(normalisation) != "character"){stop("'normalisation' must be of class 'character'")} # make sure 'normalisation' is a character
  if(!normalisation %in% c("none","sum","size class")){stop("'normalisation' must equal 'none', 'sum' or 'size class'")} # make sure normalisation equals 'none', 'sum' or 'size class'

  # clean matrix
  M[M > 0] <- 1 # ensure M is binary
  if(is.null(rownames(M))){ # if M has no row names, give it some
    rownames(M) <- paste0("r", 1:nrow(M))
  }
  if(is.null(colnames(M))){ # if M has no column names, give it some
    colnames(M) <- paste0("c", 1:ncol(M))
  }
  rn <- rownames(M) # record row names
  cn <- colnames(M) # record column names
  dimnames(M) <- NULL # strip row and column names

  # calculate inputs
  Z <- nrow(M) # number of row species
  P <- ncol(M) # number of column species

  jZ <- matrix(rep(1, Z)) # Z x 1 vector of ones
  jP <- matrix(rep(1, P)) # P x 1 vector of ones
  J <- matrix(rep(1, Z * P), nrow = Z, ncol = P) # Z x P matrix of ones
  JP <- matrix(rep(1, P * P), nrow = P, ncol = P) # P x P matrix of ones
  JZ <- matrix(rep(1, Z * Z), nrow = Z, ncol = Z) # Z x Z matrix of ones

  if(six_node == TRUE){
    JP3 <- array(rep(1, P * P * P), c(P, P, P)) # P x P x P array of ones
    JZ3 <- array(rep(1, Z * Z * Z), c(Z, Z, Z)) # Z x Z x Z array of ones
    KP3 <- JP3 # P x P x P matrix, 0 if two indices are equal, 1 otherwise
    for (i in 1 : P){
      for (j in 1 : P){
        KP3[i,j,j] <- 0
        KP3[j,i,j] <- 0
        KP3[j,j,i] <- 0
      }
    }
    KZ3 <- JZ3 # Z x Z x Z matrix, 0 if two indices are equal, 1 otherwise
    for (i in 1 : Z){
      for (j in 1 : Z){
        KZ3[i,j,j] <- 0
        KZ3[j,i,j] <- 0
        KZ3[j,j,i] <- 0
      }
    }
  }

  MT <- t(M) # transpose of M
  N <- J - M # complement of M
  NT <- t(N) # transpose of the complement of M

  dZ <- M %*% jP # degrees of the row species
  dP <- MT %*% jZ # degrees of the column species

  Z <- M %*% MT
  Y <- M %*% NT
  X <- N %*% MT

  P <- MT %*% M
  Q <- MT %*% N
  R <- NT %*% M

  if(six_node == TRUE){
    AZ <- maketensor(M, M)
    BZ <- maketensor(M, N)
    CZ <- maketensor(N, M)
    DZ <- maketensor(N, N)

    AP <- maketensor(MT, MT)
    BP <- maketensor(MT, NT)
    CP <- maketensor(NT, MT)
    DP <- maketensor(NT, NT)

    MTA <- multtensor(MT, AZ)
    MTB <- multtensor(MT, BZ)
    MTC <- multtensor(MT, CZ)
    MTD <- multtensor(MT, DZ)

    MA <- multtensor(M, AP)
    MB <- multtensor(M, BP)
    MC <- multtensor(M, CP)
    MD <- multtensor(M, DP)

    NTA <- multtensor(NT, AZ)
    NTB <- multtensor(NT, BZ)
    NTC <- multtensor(NT, CZ)

    Na <- multtensor(N, AP) # because NA already means something
    NB <- multtensor(N, BP)
    NC <- multtensor(N, CP)
  }

  # create results containers
  if(six_node == FALSE){
    pos_row <- matrix(0, ncol = 46, nrow = nrow(M), dimnames = list(rn,paste0("np",1:46)))
    pos_col <- matrix(0, ncol = 46, nrow = ncol(M), dimnames = list(cn,paste0("np",1:46)))
  } else {
    pos_row <- matrix(0, ncol = 148, nrow = nrow(M), dimnames = list(rn,paste0("np",1:148)))
    pos_col <- matrix(0, ncol = 148, nrow = ncol(M), dimnames = list(cn,paste0("np",1:148)))
  }

  # count positions
  if(six_node == FALSE){
    for(i in 1:46){
      rc <- rowcolumn(i)
      f <- countposition(M = M, p = i, jZ = jZ, jP = jP, JP = JP, JZ = JZ, MT = MT, N = N, NT = NT, dZ = dZ, dP = dP, Z = Z, Y = Y, X = X, P = P, Q = Q, R = R)
      if(rc == "row"){
        pos_row[,i] <- f
      } else {
        pos_col[,i] <- f
      }
    }
  } else {
    for(i in 1:148){
      rc <- rowcolumn(i)
      f <- countposition(M = M, p = i, jZ = jZ, jP = jP, JP = JP, JZ = JZ, JP3 = JP3, JZ3 = JZ3, KP3 = KP3, KZ3 = KZ3, MT = MT, N = N, NT = NT, dZ = dZ, dP = dP, Z = Z, Y = Y, X = X, P = P, Q = Q, R = R, MTA = MTA, MTB = MTB, MTC = MTC, MTD = MTD, MA = MA, MB = MB, MC = MC, MD = MD, NTA = NTA, NTB = NTB, NTC = NTC, Na = Na, NB = NB, NC = NC)
      if(rc == "row"){
        pos_row[,i] <- f
      } else {
        pos_col[,i] <- f
      }
    }
  }

  # normalisation
  if(normalisation != "none"){
    pos_row <- normalise_positions(pc = pos_row, type = normalisation)
    pos_col <- normalise_positions(pc = pos_col, type = normalisation)
  }

  # output
  if(level == "all"){
    out <- rbind(pos_row, pos_col)
    out <- as.data.frame(out)
    return(out)
  } else if(level == "rows"){
    pos_row <- as.data.frame(pos_row)
    return(pos_row)
  } else if(level == "columns"){
    pos_col <- as.data.frame(pos_col)
    return(pos_col)
  }
}
