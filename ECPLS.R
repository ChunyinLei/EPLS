# The main function of Elliptical Covariance Partial Least Squares (ECPLS)
ECPLS.fit <- function (X, Y, ncomp, center = TRUE, stripped = FALSE, ...) {
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
  
  for (a in 1:ncomp) {
    Kendall_xy <- Multi_Kendall_tau(X, Y)
    sigma_xy <- Robust_Sigma(X, Xmeans) %*% Kendall_xy %*% t(Robust_Sigma(Y, Ymeans))
    
    r.a <- svd(sigma_xy)$u[, 1]
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
    sigma_xy <- sigma_xy - v.a %*% crossprod(v.a, sigma_xy)
    
    # deflation
    X <- X - t.a %*% t(p.a)
    Y <- Y - t.a %*% t(q.a)
    
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
    residuals <- - fitted + c(Y)
    fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
    
    ## Add dimnames and classes:
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
