countmotif <- function(x, motif, z, p, Tz, Tp, lTz, lTp, JP = NULL, JZ = NULL, P = NULL, Q = NULL, R = NULL, Z = NULL, Y = NULL, X = NULL, dP = NULL, jP = NULL, dZ = NULL, jZ = NULL, J3 = NULL, MA = NULL, MB = NULL, MC = NULL, MD = NULL, Na = NULL, NB = NULL, NC = NULL){
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
  } else if(motif == 18){
    sum(dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) * (dP - 4 * jP)) / 120
  } else if(motif == 19){
    sum(P * Q * (Q - JP) * (Q - 2 * JP)) / 6
  } else if(motif == 20){
    sum(Q * P * R * (Q - JP)) / 2
  } else if(motif == 21){
    sum(Q * (Q - JP) * P * (P - JP)) / 4
  } else if(motif == 22){
    sum(Q * P * (P - JP) * R) / 4
  } else if(motif == 23){
    sum(P * (P - JP) * (P - 2 * JP) * Q) / 6
  } else if(motif == 24){
    sum(P * (P - JP) * (P - 2 * JP) * (P - 3 * JP)) / 48 - sum(dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 48
  } else if(motif == 25){
    sum(MA * NB * (NB - J3)) / 4
  } else if(motif == 26){
    if(z <= p){
      sum(MD * MB * MC) / 2
    } else {
      sum(NB * NC * MA) / 2
    }
  } else if(motif == 27){
    if(z <= p){
      sum(NB * NC * MA) / 2
    } else {
      sum(MD * MB * MC) / 2
    }
  } else if(motif == 28){
    sum(MB * MC * NC)
  } else if(motif == 29){
    sum(MA * MB * MD)
  } else if(motif == 30){
    if(z <= p){
      sum(MA * MB * NC) / 2
    } else {
      sum(MB * (MB - J3) * MC) / 2
    }
  } else if(motif == 31){
    if(z <= p){
      sum(MA * NB * (MA - J3)) / 4
    } else {
      sum(MB * (MB - J3) * MA) / 4
    }
  } else if(motif == 32){
    if(z <= p){
      sum(MB * (MB - J3) * MC) / 2
    } else {
      sum(MA * MB * NC) / 2
    }
  } else if(motif == 33){
    if(z <= p){
      sum(MB * (MB - J3) * MA) / 4
    } else {
      sum(MA * NB * (MA - J3)) / 4
    }
  } else if(motif == 34){
    sum(Na * MB * MC) / 6
  } else if(motif == 35){
    sum(MA * MB * MC) / 2
  } else if(motif == 36){
    sum(MA * (MA - J3) * MB) / 4
  } else if(motif == 37){
    sum(MA * (MA - J3) * (MA - 2 * J3)) / 36
  } else if(motif == 38){
    sum(Z * Y * (Y - JZ) * (Y - 2 * JZ)) / 6
  } else if(motif == 39){
    sum(Z * (Y - JZ) * X * Y) / 2
  } else if(motif == 40){
    sum(Z * (Z - JZ) * Y * (Y - JZ)) / 4
  } else if(motif == 41){
    sum(Z * (Z - JZ) * X * Y) / 4
  } else if(motif == 42){
    sum(Z * (Z - JZ) * (Z - 2 * JZ) * Y) / 6
  } else if(motif == 43){
    sum(Z * (Z - JZ) * (Z - 2 * JZ) * (Z - 3 * JZ)) / 48 - sum(dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 48
  } else if(motif == 44){
    sum(dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) * (dZ - 4 * jZ)) / 120
  }
}
