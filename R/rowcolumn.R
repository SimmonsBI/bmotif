rowcolumn <- function(position){
  if(position %in% c(2,4,6,8,11,12,14,16,18,21,22,24,25,28,29,31,34,35,38,41,42,44,46)){
    return("row")
  } else{
    return("column")
  }
}