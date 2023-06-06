#updater functions

b_update <- function(X, Y, Z, beta, bsigma, Sigma, XtX) {
  tmp_b <- matrix(nrow = ncol(X), ncol = ncol(Y))
  for(colx in 1:ncol(X)) {
    m_gj <- as.numeric(solve(XtX[colx] + 1/bsigma)) * (t(X[,colx]) %*% (Y - eigenMapMatMult(X[,-colx], beta[-colx,])))
    m_gj <- Z[colx,] * m_gj
    sig_gj <- (Z[colx,] %*% t(Z[colx,])) * (as.numeric(solve(XtX[colx] + (1/bsigma)))) * Sigma + (1 - Z[colx,] %*% t(Z[colx,])) * (bsigma * Sigma)
    sig_gj[which(Z[colx,]==1),which(Z[colx,]!=1)] <- sqrt((1/bsigma)  * as.numeric(solve(XtX[colx] + 1/bsigma))) * (bsigma * Sigma)[which(Z[colx,]==1),which(Z[colx,]!=1)]
    sig_gj[which(Z[colx,]!=1),which(Z[colx,]==1)] <- sqrt((1/bsigma)  * as.numeric(solve(XtX[colx] + 1/bsigma))) * (bsigma * Sigma)[which(Z[colx,]!=1),which(Z[colx,]==1)]
    tmp_b[colx,] <- MASS::mvrnorm(n = 1, mu = m_gj, Sigma = sig_gj)
  }
  return(tmp_b)
}

bsigma_update <- function(X, Y, b, Sigma) {
  return(LaplacesDemon::rinvgamma(n = 1, shape = 0.5 * ncol(X) * ncol(Y), scale = 0.5 * pracma::Trace(b %*% solve(Sigma) %*% t(b))))
}

gpi_update <- function(alpha) {
  return(rbeta(n = 1, shape1 = 1 + sum(alpha), shape2 = 1 + sum(1-alpha)))
}

chi_update <- function(gamma, AMAT, a) {
  mu <- c(AMAT %*% a)
  potential_chis <- rnorm(length(gamma), mu, sd = 1)
  return(ifelse(gamma == 0, -abs(potential_chis), abs(potential_chis)))
}

a_update <- function(AMAT, chi) {
  return(c(MASS::mvrnorm(n = 1, mu = solve(t(AMAT) %*% AMAT + diag(2) * (1/5)) %*% t(AMAT) %*% chi, Sigma = solve(t(AMAT) %*% AMAT + diag(2) * (1/5)))))
}

spi_update <- function(AMAT, a) {
  return(pnorm(0, mean = AMAT %*% a, sd = 1, lower.tail = FALSE))
}

tpi_update <- function(omega) {
  return(rbeta(n = 1, shape1 = 1 + sum(omega), shape2 = 1 + sum(1-omega)))
}

alpha_update <- function(X, Y, b, group, pi0, alpha, gamma, omega, Sigma) {
  gunique <- unique(group)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  for(g in 1:length(gunique)) {
    amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
    gg <- which(group == g)
    amatg <- 0 * amat
    amatg[gg,] <- 1
    amatng <- amat
    amatng[gg,] <- 0
    Zg <- amatg * gmat * omat
    Zng <- amatng * gmat * omat
    Yng <- Y - eigenMapMatMult(X, (Zng * b))
    Bg <- Zg * b
    XB <- eigenMapMatMult(X, Bg)
    calc1 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Yng - XB, SSigma), t(Yns - XB)))
    calc2 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Yng, SSigma), t(Yng)))
    prob <- pi0 / (pi0 + (1-pi0) * exp(calc2 - calc1))
    alpha[g] <- rbinom(1, 1, prob)
  }
  return(alpha)
}

gamma_update <- function(X, Y, b, group, pi1, alpha, gamma, omega, Sigma) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  for(xcol in 1:ncol(X)) {
    gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
    gmats <- gmat * 0
    gmats[xcol,] <- 1
    gmatns <- gmat
    gmatns[xcol,] <- 0
    Zs <- amat * gmats * omat
    Zns <- amat * gmatns * omat
    Yns <- Y - eigenMapMatMult(X, (Zns * b))
    Bs <- Zs * b
    XB <- eigenMapMatMult(X, Bs)
    calc1 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Yns - XB, SSigma), t(Yns - XB)))
    calc2 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Yns, SSigma), t(Yns)))
    prob <- pi1[xcol] / (pi1[xcol] + (1-pi1[xcol]) * exp(calc2 - calc1))
    gamma[xcol] <- rbinom(1, 1, prob)
  }
  return(gamma)
}

omega_update <- function(X, Y, b, group, pi2, alpha, gamma, omega, Sigma) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  tracker <- 1
  for(ycol in 1:ncol(Y)) {
    for(xcol in 1:ncol(X)) {
      omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
      omatt <-  0 * omat
      omatt[xcol,ycol] <- 1
      omatnt <- omat
      omatnt[xcol,ycol] <- 0
      Zt <- amat * gmat * omatt
      Znt <- amat * gmat * omatnt
      Ynt <- Y - eigenMapMatMult(X, (Znt * b))
      Bt <- Zt * b
      XB <- eigenMapMatMult(X, Bt)
      calc1 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Ynt - XB, SSigma), t(Ynt - XB)))
      calc2 <- (-0.5) * pracma::Trace(eigenMapMatMult(eigenMapMatMult(Ynt, SSigma), t(Ynt)))
      prob <- pi2 / (pi2 + (1-pi2) * exp(calc2 - calc1))
      omega[tracker] <- rbinom(1,1, prob)
      tracker <- tracker + 1
    }
  }
  return(omega)
}

Z_update <- function(X, Y, group, alpha, gamma, omega) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  return(amat * gmat * omat)
}

beta_update <- function(b, Z) {
  return(Z * b)
}

Sigma_update <- function(X, Y, beta, b, bsigma) {
  return(LaplacesDemon::rinvwishart(nu = ncol(X) + ncol(Y) + nrow(Y), S = diag(ncol(Y)) + t(Y - eigenMapMatMult(X, beta)) %*% (Y - eigenMapMatMult(X, beta)) + (1/bsigma) * eigenMapMatMult(t(b), b)))
}
