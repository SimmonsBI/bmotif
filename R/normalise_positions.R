normalise_positions <- function(pc,type){
  if(!ncol(pc) %in% c(46, 148)){stop("Something has gone very wrong: pc does not have 46 or 148 columns")}
  if(type == "size class"){
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:6)] <- pc[,paste0("np",3:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:6)]))
    pc[,paste0("np",7:16)] <- pc[,paste0("np",7:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:16)]))
    pc[,paste0("np",17:46)] <- pc[,paste0("np",17:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:46)]))
    if(ncol(pc) == 148){
      pc[,paste0("np",47:148)] <- pc[,paste0("np",47:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:148)]))
    }
  } else if(type == "sum"){
    pc[,paste0("np",1:ncol(pc))] <- pc[,paste0("np",1:ncol(pc))]/sapply(1:nrow(pc), function(x) sum(pc[x,]))
  } else {
    stop("'type' must be a character string equal to 'sum' or 'size class'") # if 'type' does not equal 'sum' or 'size class' return an error
  }
  return(pc)
}
