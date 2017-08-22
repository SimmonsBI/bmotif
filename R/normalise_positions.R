normalise_positions <- function(pc,type){
  if(type == "size class"){
    pc[,paste0("p",1:2)] <- pc[,paste0("p",1:2)]/apply(pc[,paste0("p",1:2)], 1, sum)
    pc[,paste0("p",3:6)] <- pc[,paste0("p",3:6)]/apply(pc[,paste0("p",3:6)], 1, sum)
    pc[,paste0("p",7:16)] <- pc[,paste0("p",7:16)]/apply(pc[,paste0("p",7:16)], 1, sum)
    pc[,paste0("p",17:46)] <- pc[,paste0("p",17:46)]/apply(pc[,paste0("p",17:46)], 1, sum)
    if(ncol(pc) == 148){
      pc[,paste0("p",47:148)] <- pc[,paste0("p",47:148)]/apply(pc[,paste0("p",47:148)], 1, sum)
    }
  } else if(type == "sum"){
    pc[,paste0("p",1:ncol(pc))] <- pc[,paste0("p",1:ncol(pc))]/apply(pc[,paste0("p",1:ncol(pc))], 1, sum)
  } else {
    stop("'type' must be a character string equal to 'sum' or 'size class'") # if 'type' does not equal 'sum' or 'size class' return an error
  }
  return(pc)
}
