countmotif <- function(x, motif, z, p, Tz, Tp, lTz, lTp){
  if(motif == 1){
    f <- sum(x)
    return(f)
  } else if(motif == 2){
    f <- 0
    for(i in 1:z){
      if(sum(x[i,]) >= 2){
        f <- f + choose(sum(x[i,]),2)
      }
    }
    return(f)
  } else if(motif == 3){
    f <- 0
    for(i in 1:p){
      if(sum(x[,i]) >= 2){
        f <- f + choose(sum(x[,i]),2)
      }
    }
    return(f)
  } else if(motif == 4){
    f <- 0
    for(i in 1:p){
      if(sum(x[,i]) >= 3){
        f <- f + choose(sum(x[,i]),3)
      }
    }
    return(f)
  } else if(motif == 5){
    f <- 0
    f <- sum(diag(x %*% t(x) %*% x %*% t(1-x)))
    return(f)
  } else if(motif == 6){
    f1 <- 0
    for(i in 1:z){
      if(i < z){
        for(j in (i+1):z){
          if(Tz[i,j] >= 2){
            f1 <- f1 + choose(Tz[i,j],2)
          }
        }
      }
    }
    f2 <- 0
    for(i in 1:p){
      if(i < p){
        for(j in (i+1):p){
          if(Tp[i,j] >= 2){
            f2 <- f2 + choose(Tp[i,j],2)
          }
        }
      }
    }
    if(f1 == f2){
      return(f1)
    } else {
      return(0)
    }
  } else if(motif == 7){
    f <- 0
    for(i in 1:z){
      if(sum(x[i,]) >= 3){
        f <- f + choose(sum(x[i,]),3)
      }
    }
    return(f)
  } else if(motif == 8){
    f <- 0
    for(i in 1:p){
      if(sum(x[,i]) >= 4){
        f <- f + choose(sum(x[,i]),4)
      }
    }
    return(f)
  } else if(motif == 9){
    f <- 0
    for(i in 1:lTp){
      for(i2 in 1:lTp){
        if(Tp[i,i]-Tp[i,i2]>=2){
          f <- f + as.numeric(choose((Tp[i,i]-Tp[i,i2]),2) %*% Tp[i,i2])
        }
      }
    }
    return(f)
  } else if(motif == 10){
    f <- 0
    for(i in 1:lTp){
      if(i < lTp){
        for(i2 in (i+1):lTp){
          f <- f + as.numeric((Tp[i,i] - Tp[i,i2]) %*% (Tp[i2,i2] - Tp[i,i2]) %*% Tp[i,i2])
        }
      }
    }
    return(f)
  } else if(motif == 11){
    f <- 0
    for(i in 1:lTp){
      for(i2 in 1:lTp){
        if(Tp[i,i2]>=2){
          f <- f + as.numeric(choose((Tp[i,i2]),2) %*% (Tp[i,i]-Tp[i,i2]))
        }
      }
    }
    return(f)
  } else if(motif == 12){
    f <- 0
    for(i in 1:lTp){
      if(i<lTp){
        for(j in (i+1):lTp){
          if(Tp[i,j]>=3){
            f <- f + choose(Tp[i,j],3)
          }
        }
      }
    }
    return(f)
  } else if(motif == 13){
    f <- 0
    for(i in 1:lTz){
      for(i2 in 1:lTz){
        if(Tz[i,i]-Tz[i,i2]>=2){
          f <- f + as.numeric(choose((Tz[i,i]-Tz[i,i2]),2) %*% Tz[i,i2])
        }
      }
    }
    return(f)
  } else if(motif == 14){
    f <- 0
    for(i in 1:lTz){
      if(i < lTz){
        for(i2 in (i+1):lTz){
          f <- f + as.numeric((Tz[i,i] - Tz[i,i2]) %*% (Tz[i2,i2] - Tz[i,i2]) %*% Tz[i,i2])
        }
      }
    }
    return(f)
  } else if(motif == 15){
    f <- 0
    for(i in 1:lTz){
      for(i2 in 1:lTz){
        if(Tz[i,i2]>=2){
          f <- f + as.numeric(choose((Tz[i,i2]),2) %*% (Tz[i,i]-Tz[i,i2]))
        }
      }
    }
    return(f)
  } else if(motif == 16){
    f <- 0
    for(i in 1:lTz){
      if(i<lTz){
        for(j in (i+1):lTz){
          if(Tz[i,j]>=3){
            f <- f + choose(Tz[i,j],3)
          }
        }
      }
    }
    return(f)
  } else if(motif == 17){
    f <- 0
    for(i in 1:z){
      if(sum(x[i,]) >= 4){
        f <- f + choose(sum(x[i,]),4)
      }
    }
    return(f)
  }
}
