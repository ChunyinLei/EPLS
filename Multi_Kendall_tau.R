## Multivariate Kendall’s tau
Multi_Kendall_tau <- function(X, Y){
  if(!is.matrix(X) && !is.data.frame(X)){
    stop("X must be either numeric vector, matrix or data.frame.")}
  if(!is.matrix(Y) && !is.data.frame(Y)){
    stop("Y must be either numeric vector, matrix or data.frame.")}
  
  n <- nrow(X)
  s <- 0
  for(j in 1:(n - 1)){
    for(i in (j + 1):n){
      if(norm((X[i, ] - X[j, ]), '2') == 0 || norm((Y[i, ] - Y[j, ]), '2') == 0){
        s <- s + matrix(0, nrow = ncol(X), ncol = ncol(Y))
      }else{
        s <- s + (as.matrix(X[i, ] - X[j, ]) %*% t(Y[i, ] - Y[j, ])) / (norm((X[i, ] - X[j, ]), '2') * norm((Y[i, ] - Y[j, ]), '2'))
      }
    }
  }
  ret <- s / choose(n, 2)
  
  return(ret)
}
