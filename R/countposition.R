countposition <- function(M, p, jZ, jP, J, JP, JZ, MT, N, NT, dZ, dP, Z, Y, X, P, Q, R){
  if (p == 1) {a <- dP}
  else if (p == 2) {a <- dZ}
  else if (p == 3) {a <- P %*% jP - dP}
  else if (p == 4) {a <- dZ * (dZ - jZ) / 2}
  else if (p == 5) {a <- dP * (dP - jP) / 2}
  else if (p == 6) {a <- Z %*% jZ - dZ}
  else if (p == 7) {a <- dP * (dP - jP) * (dP - 2 * jP) / 6}
  else if (p == 8) {a <- M %*% ((dP - jP) * (dP - 2 * jP)) / 2}
  else if (p == 9) {a <- (P * R) %*% jP}
  else if (p == 10) {a <- (P * Q) %*% jP}
  else if (p == 11) {a <- (X * Z) %*% jZ}
  else if (p == 12) {a <- (Y * Z) %*% jZ}
  else if (p == 13) {a <- (P * (P - JP)) %*% jP / 2 - dP * (dP - jP) / 2}
  else if (p == 14) {a <- (Z * (Z - JZ)) %*% jZ / 2 - dZ * (dZ - jZ) / 2}
  else if (p == 15) {a <- MT %*% ((dZ - jZ) * (dZ - 2 * jZ)) / 2}
  else if (p == 16) {a <- dZ * (dZ - jZ) * (dZ - 2 * jZ) / 6}
  else if (p == 17) {a <- dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) / 24}
  else if (p == 18) {a <- M %*% ((dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 6}
  else if (p == 19) {a <- (P * R * (R - JP)) %*% jP / 2}
  else if (p == 20) {a <- (P * Q * (Q - JP)) %*% jP / 2}
  else if (p == 21) {a <- (N * (M %*% ((Q - JP) * P))) %*% jP}
  else if (p == 22) {a <- (M * (M %*% (Q * (Q - JP)))) %*% jP / 2}
  else if (p == 23) {a <- (P * Q * R) %*% jP}
  else if (p == 24) {a <- (N * (M %*% (P * R))) %*% jP}
  else if (p == 25) {a <- (M * (M %*% (Q * R))) %*% jP / 2}
  else if (p == 26) {a <- (P * (P - JP) * R) %*% jP / 2}
  else if (p == 27) {a <- (P * (P - JP) * Q) %*% jP / 2}
  else if (p == 28) {a <- (N * (M %*% (P * (P - JP)))) %*% jP / 2}
  else if (p == 29) {a <- (M * (M %*% (Q * (P - JP)))) %*% jP}
  else if (p == 30) {a <- (P * (P - JP) * (P - 2 * JP)) %*% jP / 6 - dP * (dP - jP) * (dP - 2 * jP) / 6}
  else if (p == 31) {a <- (M * M %*% ((P - JP) * (P - 2 * JP))) %*% jP / 4 - M %*% ((dP - jP) * (dP - 2 * jP)) / 4}
  else if (p == 32) {a <- (NT * (MT %*% ((Y - JZ) * Z))) %*% jZ}
  else if (p == 33) {a <- (MT * (MT %*% (Y * (Y - JZ)))) %*% jZ / 2}
  else if (p == 34) {a <- (Z * X * (X - JZ)) %*% jZ / 2}
  else if (p == 35) {a <- (Z * Y * (Y - JZ)) %*% jZ / 2}
  else if (p == 36) {a <- (NT * (MT %*% (Z*X))) %*% jZ}
  else if (p == 37) {a <- (MT * (MT %*% (Y * X))) %*% jZ / 2}
  else if (p == 38) {a <- (Z * Y * X) %*% jZ}
  else if (p == 39) {a <- (NT * (MT %*% (Z * (Z - JZ)))) %*% jZ / 2}
  else if (p == 40) {a <- (MT * (MT %*% (Y * (Z - JZ)))) %*% jZ}
  else if (p == 41) {a <- (Z * (Z - JZ) * X) %*% jZ / 2}
  else if (p == 42) {a <- (Z * (Z - JZ) * Y) %*% jZ / 2}
  else if (p == 43) {a <- (MT * (MT %*% ((Z - JZ) * (Z - 2 * JZ)))) %*% jZ / 4 - MT %*% ((dZ - jZ) * (dZ - 2 * jZ)) / 4}
  else if (p == 44) {a <- (Z * (Z - JZ) * (Z - 2 * JZ)) %*% jZ / 6 - dZ * (dZ - jZ) * (dZ - 2 * jZ) / 6}
  else if (p == 45) {a <- MT %*% ((dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 6}
  else if (p == 46) {a <- dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) / 24}
  else {a <- -1}
  a
}