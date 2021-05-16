cp <- function(x, y, w = 1){
  
  if (length(x) != length(y)) stop('[cp] x and y must have the same length.')
  
  if (length(w) != 1 && length(w) != length(x)) stop('[cp] x and w must have the same length.')
  
  cp_x <- sum(w * x) / sum(w)
  cp_y <- sum(w * y) / sum(w)
  cp <- cbind(cp_x, cp_y)
  names(cp) <- c('cp_x', 'cp_y')
  return(cp)
}
