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
    # pc[,paste0("np",1:ncol(pc))] <- pc[,paste0("np",1:ncol(pc))]/sapply(1:nrow(pc), function(x) sum(pc[x,]))
    pc <- pc/apply(pc, 1, sum)
  } else if(type == "position"){
    pc <- apply(pc, MARGIN = 2, FUN = function(x) x/sum(x))
  } else if(type == "levelsize"){
    pc[,paste0("np",1)] <- pc[,paste0("np",1)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1)]))
    pc[,paste0("np",2)] <- pc[,paste0("np",2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",2)]))
    pc[,paste0("np",3)] <- pc[,paste0("np",3)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3)]))
    pc[,paste0("np",4)] <- pc[,paste0("np",4)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",4)]))
    pc[,paste0("np",5:6)] <- pc[,paste0("np",5:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",5:6)]))
    pc[,paste0("np",7)] <- pc[,paste0("np",7)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7)]))
    pc[,paste0("np",8)] <- pc[,paste0("np",8)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",8)]))
    pc[,paste0("np",9:12)] <- pc[,paste0("np",9:12)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",9:12)]))
    pc[,paste0("np",13:16)] <- pc[,paste0("np",13:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",13:16)]))
    pc[,paste0("np",17)] <- pc[,paste0("np",17)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17)]))
    if(ncol(pc) == 148){
      pc[,paste0("np",18)] <- pc[,paste0("np",18)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",18)]))
      pc[,paste0("np",19:24)] <- pc[,paste0("np",19:24)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",19:24)]))
      pc[,paste0("np",25:37)] <- pc[,paste0("np",25:37)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",25:37)]))
      pc[,paste0("np",38:43)] <- pc[,paste0("np",38:43)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",38:43)]))
      pc[,paste0("np",44)] <- pc[,paste0("np",44)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",44)]))
    }
  } else {
    stop("'type' must be a character string equal to 'sum', 'size class', 'position' or 'levelsize'") # if 'type' does not equal 'sum' or 'size class' return an error
  }
  return(pc)
}
