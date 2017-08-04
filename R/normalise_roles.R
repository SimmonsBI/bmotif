normalise_roles <- function(M,type){
  if(type == "within"){
    M[,paste0("p",1:2)] <- M[,paste0("p",1:2)]/apply(M[,paste0("p",1:2)], 1, sum)
    M[,paste0("p",3:6)] <- M[,paste0("p",3:6)]/apply(M[,paste0("p",3:6)], 1, sum)
    M[,paste0("p",7:16)] <- M[,paste0("p",7:16)]/apply(M[,paste0("p",7:16)], 1, sum)
    M[,paste0("p",17:46)] <- M[,paste0("p",17:46)]/apply(M[,paste0("p",17:46)], 1, sum)
  } else if(type == "across"){
    M[,paste0("p",1:46)] <- M[,paste0("p",1:46)]/apply(M[,paste0("p",1:46)], 1, sum)
  } else {
    stop("'type' must be a character string equal to 'within' or 'across'") # if 'type' does not equal 'within' or 'across' return an error
  }
  return(M)
}