context("check_motif")
### Testing if the check_motif function works the way we want it to work #####

# define functions
gen_all_mot_perm <- function(all_mot = -1) {
  # all_mot should be a list of all motifs that we want to permute
  # the function generates all permutations of the motifs
  # returns a nested list
  # all_mot_perm[[i]] is a list of all permutations of the element i in all_mot
  # if called without any argument, will load in data

  if (all_mot == -1) {
    # load('all_motifs.RData')
    all_mot <- motifs
  }

  # library(gtools)
  all_mot_perm <- list()
  v <- list()
  v[[1]] <- NA
  v[[2]] <- 1
  v[[3]] <- 2:3
  v[[4]] <- 4:7
  v[[5]] <- 8:17
  v[[6]] <- 18:44

  # loop through possible motif sizes
  for (n in 2:6) {
    perm <- list()

    # loop through poss number of vertices in first class
    for (i in 1:(n-1)) {
      # now have i vertices in first class, n-i in second
      pr <- gtools::permutations(i, i, 1:i, repeats = FALSE)
      pc <- gtools::permutations(n-i, n-i, 1:(n-i), repeats = FALSE)
      perm[[i]] <- list(pr, pc)
    }

    # now permute all the motifs
    for (j in v[[n]]) {
      mot <- all_mot[[j]]
      k <- nrow(mot)
      pr <- perm[[k]][[1]]
      pc <- perm[[k]][[2]]
      all_mot_perm[[j]] <- bipgraph_isom_class(mot, pr, pc)
    }
  }
  all_mot_perm
}

bipgraph_isom_class <- function(M, poss_perm_rows, poss_perm_cols ) {
  # M is adjacency matrix of a bipartite graph (i.e. not nec. square matrix)
  # now we want to find all bipartite graphs obtained by just permututing rows and columns
  # poss_perm_rows is a matrix with possible permutations of the row indices
  # poss_perm_cols ... same for column indices
  # we return a list with all isomorphic bipartite graphs

  # poss_perm_rows <- permutations(N, N, 1:N, repeats = FALSE)
  # poss_perm_cols <- permutations(M, M, 1:M, repeats = FALSE)

  # these are given as arguments to this function so we don't have to compute them again, every time the function is called

  I <- nrow(poss_perm_rows)
  J <- nrow(poss_perm_cols)

  res <- list()

  for (i in 1:I) {
    M_pr <- M[poss_perm_rows[i,], , drop = FALSE]
    for (j in 1:J) {
      M_pc <- M_pr[,poss_perm_cols[j,], drop = FALSE]
      len <- length(res)
      res[[len + 1]] <- M_pc
    }
  }

  res <- unique(res)
  return(res)

}

check_motif <- function(minor, all_mot_perm) {
  # all_mot_perm[[i]] is a list of all permutations of motif i
  # minor is a subgraph
  # function returns the motif type of the given minor

  # delete row and column names
  dimnames(minor) <- NULL

  # don't have to try all motifs, only the ones of same dimension
  poss_mot <- motif_dimensions(nrow(minor), ncol(minor))

  for (i in poss_mot) {
    iso_class <- all_mot_perm[[i]]
    for (j in 1:length(iso_class)) {
      candidate <- iso_class[[j]]
      dimnames(candidate) <- NULL
      if (identical(minor, candidate)) {
        return(i)
      }
    }
  }
  return(-1)
}

motif_dimensions <- function(m, n) {
  # This function takes as arguments the dimensions of a subgraph
  # It returns a vector with all possible motif types of the same dimension
  if (m == 1 & n == 1) {
    return(1)
  }
  if (m == 1 & n == 2) {
    return(2)
  }
  if (m == 2 & n == 1) {
    return(3)
  }
  if (m == 3 & n == 1) {
    return(4)
  }
  if (m == 2 & n == 2) {
    return(5:6)
  }
  if (m == 1 & n == 3) {
    return(7)
  }
  if (m == 4 & n == 1) {
    return(8)
  }
  if (m == 3 & n == 2) {
    return(9:12)
  }
  if (m == 2 & n == 3) {
    return(13:16)
  }
  if (m == 1 & n == 4) {
    return(17)
  }
  if (m == 5 & n == 1) {
    return(18)
  }
  if (m == 4 & n == 2) {
    return(19:24)
  }
  if (m == 3 & n == 3) {
    return(25:37)
  }
  if (m == 2 & n == 4) {
    return(38:43)
  }
  if (m == 1 & n == 5) {
    return(44)
  }

}

# tests
all_mot_perm <- gen_all_mot_perm()

for(i in 1:17) {
  mot <- motifs[[i]]
  N <- nrow(mot)
  M <- ncol(mot)
  pr <- gtools::permutations(N, N, 1:N, repeats = FALSE)
  pc <- gtools::permutations(M, M, 1:M, repeats = FALSE)
  # pick a random permutation for rows and cols
  rr <- floor(stats::runif(1, min=1, max= nrow(pr) + 1))
  rc <- floor(stats::runif(1, min=1, max= nrow(pc) + 1))

  mot2 <- mot[pr[rr, ], ,drop = FALSE]
  mot3 <- mot2[, pc[rc, ], drop = FALSE]

  numb <- check_motif(mot3, all_mot_perm)

  expect(numb == i)
}
