## Robust Sigma
Robust_Sigma <- function(X, Xmeans, sigma0){
  sigma_variables <- c()
  for(i in 1:ncol(X)){
    sigma_variables[i] <- var(X[, i])
  }
  
  sigma0_upper <- min(sigma_variables)
  sigma0 <- sigma0_upper/2 # sigma0 can be chosen to be any positive number less than sigma0_upper
  
  rsigma <- c()
  for(j in 1:ncol(X)){
    rsigma[j] <- max(c(sigma0, var(X[, j]) - Xmeans[j]^2))
  }
  
  robust_sigma <- diag(rsigma, ncol = ncol(X))
  return(robust_sigma)
} 
