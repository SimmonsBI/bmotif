mprofile <- function(M, normalise){
  #' Count motifs in a network
  #'
  #' This function counts the number of times each of the 17 motifs up to five nodes occurs in a network
  #' @param M A numeric matrix representing interactions between two groups of species. Each row corresponds to a species in one level
  #' and each column corresponds to a species in the other level. Elements of M are positive numbers if species do interact, and 0 
  #' otherwise. Formally, M is an incidence matrix. When species i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
  #' If interactions are weighted (matrix elements are greater than 1), the function will automatically convert the matrix to a binary
  #' matrix.
  #' @param normalise Logical; should motif frequencies be normalised to control for network size?
  #' @export
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
    out <- data.frame(motif = 1:17, nodes = c(1,rep(2,2),rep(4,4),rep(5,10)), frequency = NA, normalise_sum = NA, normalise_sizeclass = NA, normalise_speciessets = NA)
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
    
    # calculate normalised frequency as proportion of possible species sets
    sets <- species_sets(M)
    out$normalise_speciessets <- out$frequency/sets
  }
  
  # output
  return(out)
}


# mprofile <- function(network){
#   prof <- setNames(object = rep(NA, 16), nm = paste0("m",1:16))
#   for(i in 1:16){
#     prof[i] <- countmotif(network, i)
#   }
#   # sapply(1:16, function(i) countmotif(network, i))
#   return(prof)
# }

# mprofile <- function(M){
#   # clean matrix
#   M[M > 0] <- 1 # ensure M is binary
#   dimnames(M) <- NULL # strip row and column names
#   
#   # calculate inputs
#   z <- dim(M)[1]
#   p <- dim(M)[2]
#   Tz <- M %*% t(M)
#   Tp <- t(M) %*% M
#   lTz <- dim(Tz)[2]
#   lTp <- dim(Tp)[2]
#   
#   # create results container
#   prof <- setNames(object = rep(NA, 16), nm = paste0("m",1:16))
#   
#   # count motifs
#   for(i in 1:16){
#     prof[i] <- countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp)
#   }
#   # sapply(1:16, function(i) countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp))
#   
#   # assemble output
#   # out <- matrix(nrow = 17, ncol = 5, dimnames = list(NULL, c("motif","nodes","frequency","norm_freq_sizeclass","norm_freq_speciessets")))
#   out <- data.frame(motif = 1:17, nodes = c(1,rep(2,2),rep(4,4),rep(5,10)), frequency = NA, norm_freq_sizeclass = NA, norm_freq_speciessets = NA)
#   out$frequency <- sapply(1:17, function(i) countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp))
#   
#   return(prof)
# }

# microbenchmark(
#   for(i in 1:16){
#     prof[i] <- countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp)
#   }
# , times = 100)
# 
# microbenchmark(
#   out[,"frequency"] <- sapply(1:17, function(i) countmotif(x = M, motif =  i, z = z, p = p, Tz = Tz, Tp = Tp, lTz = lTz, lTp = lTp))
# , times = 100)
