## Multivariate Kendall’s tau
Multi_Kendall_tau <- function(x){
  if(!is.matrix(x) && !is.data.frame(x)){
    stop("x must be either numeric vector, matrix or data.frame.")}
  n <- nrow(x)
  s <- 0
  for (j in 1:(n - 1)){
    for (i in (j + 1):n){
      s <- s + ((x[i, ] - x[j, ]) %*% t(x[i, ] - x[j, ])) / (norm((x[i, ] - x[j, ]), '2') ^ 2)
      ret <- s / choose(n, 2)
    }
  }
  return(ret)
}
