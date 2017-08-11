positions <- function(M, six_node, level = "all", normalisation = "none"){
  #' Calculate node position vectors
  #'
  #' Counts the number of times each node in a network occurs in each of the 46 positions found within the 17 motifs up to five nodes
  #' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
  #' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes do interact, and 0
  #' otherwise. Formally, M is an incidence matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
  #' If interactions are weighted (non-zero matrix elements can take values greater than 1), the function will automatically convert the matrix to a binary
  #' matrix.
  #' @param six_node Logical; should six node motifs be counted?
  #' @param  level Which node level should positions be calculated for: \code{rows}, \code{columns} or \code{all}?  Defaults to \code{all}.
  #' @param normalisation Which normalisation should be used: \code{none}, \code{across} or \code{within}?  Defaults to \code{none}.
  #' @details The \code{level} argument controls which node group positions are calculated for. \code{rows} returns position counts for all nodes in rows, \code{columns}
  #' returns counts for all nodes in columns, and \code{all} return counts for all nodes in the network.
  #'
  #' Nodes with more interactions will tend to appear in more positions. Normalisation helps control for this.
  #' \code{none} performs no normalisation and will return the raw position counts.
  #' \code{across} divides position counts for each node by the total number of times that node appears in any position.
  #' \code{within} divides position counts for each node by the total number of times that node appears in any position within the same motif size class.
  #'
  #' Warning: including six node motifs is fine for most networks. However, for large networks, counting six node motifs can be slow and memory intensive. In some cases, R can crash if there is not enough memory.
  #' @return
  #' By default, returns a data frame with 46 columns, one for each motif position. If \code{six_node} is TRUE, there will be 144 columns (one for each position).
  #' For a network with A rows and P columns, by default (where \code{level} = "all") the data frame has A + P rows, one for each node. If \code{level} = "rows", the data frame will have A rows, one for each row node;
  #' if \code{level} = "columns", it will have P rows, one for each column node.
  #'
  #' By default, the elements of this data frame will be the raw position counts. If \code{normalisation} is set to "within" or "across", the elements will be
  #' the normalised position counts.
  #' @export
  #' @examples
  #' set.seed(123)
  #' row <- 15
  #' col <- 15
  #' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
  #' positions(M = m, six_node = TRUE, level = "all", normalisation = "none")

  # check inputs
  if(class(M) != "matrix"){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length 0")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(class(level) != "character"){stop("'level' must be of class 'character'")} # make sure level is a character
  if(!level %in% c("rows","columns","all")){stop("'level' must equal 'rows', 'columns' or 'all'")} # make sure level equals 'rows', 'columns' or 'all'
  if(class(normalisation) != "character"){stop("'normalisation' must be of class 'character'")} # make sure 'normalisation' is a character
  if(!normalisation %in% c("none","within","across")){stop("'normalisation' must equal 'none', 'within' or 'across'")} # make sure normalisation equals 'none', 'within' or 'across'

  # clean matrix
  M[M > 0] <- 1 # ensure M is binary
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
    pos_row <- matrix(0, ncol = 46, nrow = nrow(M), dimnames = list(rn,paste0("p",1:46)))
    pos_col <- matrix(0, ncol = 46, nrow = ncol(M), dimnames = list(cn,paste0("p",1:46)))
  } else {
    pos_row <- matrix(0, ncol = 148, nrow = nrow(M), dimnames = list(rn,paste0("p",1:148)))
    pos_col <- matrix(0, ncol = 148, nrow = ncol(M), dimnames = list(cn,paste0("p",1:148)))
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
    pos_row <- normalise_positions(M = pos_row, type = normalisation)
    pos_col <- normalise_positions(M = pos_col, type = normalisation)
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
