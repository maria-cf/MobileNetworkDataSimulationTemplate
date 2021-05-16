rd <- function(x, y, w = 1, center = NULL, method = 'euclidean', p = 2){
  
  n <- length(x)
  
  if (n != length(y)) stop('[rd] x and y must have the same length.')
  
  if (length(w) != 1 && length(w) != n) stop('[rd] x and w must have the same length.')
  
  if (is.null(center)) {
    
    center <- cp(x, y, w) 
    names(center) <- c('cp_x', 'cp_y')
  
  }

  coordMatrix <- rbind(center, cbind(x, y))
  distance <- dist(coordMatrix, method = method, p = p)[1:n]
  w <- w / sum(w)
  rd <- sqrt(sum(w * distance**2))
  return(rd)  
  
}
