library(pls)
library(rpls)
library(sgt)
library(actuar)
library(nipals)
library(LaplacesDemon)
library(NonNorMvtDist)

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

Multi_kenpls.fit0 <- function (X, Y, ncomp, center = TRUE, stripped = FALSE, ...) {
  Y <- as.matrix(Y)
  if (!stripped) {
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }
  dimnames(X) <- dimnames(Y) <- NULL
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  V <- R <- matrix(0, nrow = npred, ncol = ncomp)
  tQ <- matrix(0, nrow = ncomp, ncol = nresp)
  B <- array(0, dim = c(npred, nresp, ncomp))
  if (!stripped) {
    P <- R
    U <- TT <- matrix(0, nrow = nobj, ncol = ncomp)
    fitted <- array(0, dim = c(nobj, nresp, ncomp))
  }
  if (center) {
    Xmeans <- colMeans(X)
    X <- X - rep(Xmeans, each = nobj)
    Ymeans <- colMeans(Y)
    Y <- Y - rep(Ymeans, each = nobj)
  }else {
    Xmeans <- rep_len(0, npred)
    Ymeans <- rep_len(0, nresp)
  }
  
  XY <- cbind(X, Y)
  sigmaXY <- Multi_Kendall_tau(XY)
  m <- npred+1
  o <- npred+nresp
  sigma_xx <- sigmaXY[1:npred,1:npred]
  sigma_yy <- sigmaXY[m:o, m:o]
  sigma_xy <- sigmaXY[1:npred, m:o]
  S <- sigma_xy
  
  for (a in 1:ncomp) {
    if (nresp == 1) {
      q.a <- matrix(c(1),1,1)
    }else {
      if (nresp < npred) {
        q.a <- eigen(crossprod(S), symmetric = TRUE)$vectors[, 1]
      }else {
        q.a <- c(crossprod(S, eigen(S %*% t(S), symmetric = TRUE)$vectors[, 1]))
        q.a <- q.a/sqrt(c(crossprod(q.a)))
      }
    }
    r.a <- S %*% q.a
    t.a <- X %*% r.a
    if (center) {t.a <- t.a - mean(t.a)}
    tnorm <- sqrt(c(crossprod(t.a)))
    t.a <- t.a/tnorm
    r.a <- r.a/tnorm
    p.a <- crossprod(X, t.a)
    q.a <- crossprod(Y, t.a)
    v.a <- p.a
    if (a > 1) {v.a <- v.a - V %*% crossprod(V, p.a)}
    v.a <- v.a/sqrt(c(crossprod(v.a)))
    S <- S - v.a %*% crossprod(v.a, S)
    R[, a] <- r.a
    tQ[a, ] <- q.a
    V[, a] <- v.a
    B[, , a] <- R[, 1:a, drop = FALSE] %*% tQ[1:a, , drop = FALSE]
    if (!stripped) {
      u.a <- Y %*% q.a
      if (a > 1) 
        u.a <- u.a - TT %*% crossprod(TT, u.a)
      P[, a] <- p.a
      TT[, a] <- t.a
      U[, a] <- u.a
      fitted[, , a] <- TT[, 1:a] %*% tQ[1:a, , drop = FALSE]
    }
  }
  if (stripped) {
    list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
  }else {
    residuals <- -fitted + c(Y)
    fitted <- fitted + rep(Ymeans, each = nobj)
    objnames <- dnX[[1]]
    if (is.null(objnames)) 
      objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
    dimnames(R) <- dimnames(P) <- list(prednames, compnames)
    dimnames(tQ) <- list(compnames, respnames)
    dimnames(B) <- list(prednames, respnames, nCompnames)
    dimnames(fitted) <- dimnames(residuals) <- list(objnames, respnames, nCompnames)
    class(TT) <- class(U) <- "scores"
    class(P) <- class(tQ) <- "loadings"
    
    list(coefficients = B, scores = TT, loadings = P, Yscores = U, 
         Yloadings = t(tQ), projection = R, Xmeans = Xmeans, 
         Ymeans = Ymeans, fitted.values = fitted, residuals = residuals, 
         Xvar = colSums(P * P), Xtotvar = sum(X * X))
  }
}

method1 <- c()
method2 <- c()
method3 <- c()
method4 <- c()
pre_method1 <- c()
pre_method2 <- c()
pre_method3 <- c()
pre_method4 <- c()

t1 <- Sys.time()

for(t in 1:500){
  set.seed(t)
  
  n <- 400 # Sample size
  k <- 30 # Number of independent variables
  A <- 4 # Number of components
  # SNR <- 10 # Signal noise ratio

  # mu <- rep(0, k)
  # Sigma <- diag(k)
  # X <- matrix(rmvn(n, mu, Sigma),n,k) # Multivariate Norm

  # mu <- rep(0, k)
  # Sigma <- diag(k)
  # X <- matrix(rmvl(n, mu, Sigma),n,k) # Multivariate Laplace

  # X <- matrix(rmvlogis(n, parm1 = rep(0.2, k), parm2 = rep(0.2, k)),n,k) # Multivariate Logistic
  
  mu <- rep(0, k)
  Sigma <- diag(k)
  X <- matrix(rmvt(n, mu, Sigma, df=2),n,k) # Multivariate t
  
  # mu <- rep(0, k)
  # Sigma <- diag(k)
  # X <- matrix(rmvc(n, mu, Sigma),n,k) # Multivariate cauchy
  
  beta <- c(rnorm(k,0,0.1))
  
  e0 <- c(rnorm(n,0,1))
  # e0 <- c(rnorm(n,0,sqrt(var(X %*% beta) / SNR)))
  e1 <- c(rnorm(ceiling(0.9*n),0,1),rnorm(floor(0.1*n),10,1))
  e2 <- c(rnorm(ceiling(0.9*n),0,1),rnorm(floor(0.1*n),6,2))
  e3 <- c(rlaplace(n, location = 0, scale = 1)) # Laplace
  e4 <- c(rt(n, 3)) # t3
  e5 <- c(rt(n, 2)) # t2
  e6 <- c(rmvlogis(n, parm1 = 0.2, parm2 = 0.2)) # Logistic
  e7 <- c(rcauchy(n, location = 0, scale = 1)) # cauchy
  
  y <- X %*% beta + e5
  d <- as.data.frame(X)
  d <- cbind(y,d)
  
  #### kendall pls
  l1 <- Multi_kenpls.fit0(X,y,A)
  l1$coefficients
  
  #### simpls
  pls2 <- simpls.fit(X,y,A)
  pls2$coefficients 
  
  #### PRM
  res <- PRM(y~., data=d, a=A, wfunX=list("Fair",0.95),
             wfunY=list("Fair",0.95), center.type = "median",
             scale.type = "no",usesvd = FALSE, numit = 100, prec = 0.01)
  
  #### nipls
  l2 <- nipals(X, maxiter=1000)
  l2_coefficients <- l2$loadings[,1:A] %*% l2$R2[1:A]
  
  #### beta's mse
  method1[t] <- sum((l1$coefficients[,,A]-beta)^2/length(l1$coefficients[,,A]))
  method2[t] <- sum((res$coef-beta)^2/length(res$coef))
  method3[t] <- sum((pls2$coefficients[,,A]-beta)^2/length(pls2$coefficients[,,A]))
  method4[t] <- sum((l2_coefficients-beta)^2/length(l2_coefficients))
  
  ### pre_y's mse
  pre_method1[t] <- mean((X %*% l1$coefficients[,,A] - y)^2)
  pre_method2[t] <- mean((X %*% res$coef - y)^2)
  pre_method3[t] <- mean((X %*% pls2$coefficients[,,A] - y)^2)
  pre_method4[t] <- mean((X %*% l2_coefficients - y)^2)
}

t2 <- Sys.time()
t2 - t1 # Time difference

mean_mse_method1 <- round(mean(method1), 4)
std_mse_method1 <- round(sqrt(var(method1)), 4)
mean_mse_method2 <- round(mean(method2), 4)
std_mse_method2 <- round(sqrt(var(method2)), 4)
mean_mse_method3 <- round(mean(method3), 4)
std_mse_method3 <- round(sqrt(var(method3)), 4)
mean_mse_method4 <- round(mean(method4), 4)
std_mse_method4 <- round(sqrt(var(method4)), 4)

mean_mse_pre_method1 <- round(mean(pre_method1), 4)
std_mse_pre_method1 <- round(sqrt(var(pre_method1)), 4)
mean_mse_pre_method2 <- round(mean(pre_method2), 4)
std_mse_pre_method2 <- round(sqrt(var(pre_method2)), 4)
mean_mse_pre_method3 <- round(mean(pre_method3), 4)
std_mse_pre_method3 <- round(sqrt(var(pre_method3)), 4)
mean_mse_pre_method4 <- round(mean(pre_method4), 4)
std_mse_pre_method4 <- round(sqrt(var(pre_method4)), 4)

## result
matrix(c(mean_mse_method1,mean_mse_method2,mean_mse_method3,mean_mse_method4,
         std_mse_method1,std_mse_method2,std_mse_method3,std_mse_method4), nrow = 4, 
       dimnames = list(c("kendall pls","PRM","simpls",'nipls'),c("mean(beta's mse)","std(beta's mse)")))
       
matrix(c(mean_mse_pre_method1,mean_mse_pre_method2,mean_mse_pre_method3,mean_mse_pre_method4,
         std_mse_pre_method1,std_mse_pre_method2,std_mse_pre_method3,std_mse_pre_method4), nrow = 4, 
       dimnames = list(c("kendall pls","PRM","simpls",'nipls'),c("mean(pre_y's mse)","std(pre_y's mse)")))