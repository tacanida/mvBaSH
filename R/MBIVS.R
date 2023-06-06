#This is the primary function
#' @import Rcpp

sourceCpp("/src/matfuns.cpp")

f <- function(X, Y, group, A, burnin = 2500, niter = 5000) {
  pi0 <- numeric(burnin+niter); pi0[1] <- rbeta(1,1,1)
  pi1 <- matrix(NA, nrow = burnin+niter, ncol = length(group)); pi1[1,] <- rbeta(n = length(group), 1, 1)
  pi2 <- numeric(burnin+niter); pi2[1] <- rbeta(1,1,1)
  alpha <- matrix(NA, nrow = burnin+niter, ncol = length(unique(group))); alpha[1,] <- rbinom(n = length(unique(group)), size = 1, prob = pi0[1])
  gamma <- matrix(NA, nrow = burnin+niter, ncol = length(group)); gamma[1,] <- rbinom(n = length(group), size = 1, prob = pi1[1,])
  omega <- matrix(NA, nrow = burnin+niter, ncol = ncol(X) * ncol(Y)); omega[1,] <- rbinom(n = ncol(X) * ncol(Y), size = 1, prob = pi2[1])
  chi <- matrix(NA, nrow = burnin+niter, ncol = length(group)); chi[1,] <- rnorm(n = length(group))
  a <- matrix(NA, nrow = burnin + niter, ncol = 2); a[1,] <- rnorm(n = 2)
  Z <- array(NA, dim = c(ncol(X), ncol(Y), burnin+niter))
  Z[,,1] <- 0
  Z[which(A==1),,1] <- 1
  bsigma <- numeric(burnin+niter); bsigma[1] <- runif(1, 0, 1)
  b <- array(NA, dim = c(ncol(X), ncol(Y), burnin+niter))
  b[,,1] <- solve(t(X) %*% X) %*% t(X) %*% Y #least squares starting point
  beta <-  array(NA, dim = c(ncol(X), ncol(Y), burnin+niter))
  beta[,,1] <- Z[,,1] * b[,,1]
  Sigma <- array(NA, dim = c(ncol(Y), ncol(Y), burnin+niter))
  Sigma[,,1] <- t(Y - X %*% beta[,,1]) %*% (Y - X %*% beta[,,1]) / nrow(Y-1)
  n <- nrow(X)
  p <- ncol(X)
  X <- matrix(as.double(X), nrow = n, ncol = p, byrow = FALSE)
  AMAT <- cbind(1, A)
  XtXgj <- apply(X, 2, FUN = function(x) sum(x^2))

  for(i in 2:(burnin+niter)) {
    b[,,i] <- b_update(X = X, Y = Y, Z = Z[,,i-1], beta = beta[,,i-1], bsigma = bsigma[i-1], Sigma = Sigma[,,i-1], XtX = XtXgj)
    bsigma[i] <- bsigma_update(X = X, Y = Y, b = b[,,i], Sigma = Sigma[,,i-1])
    pi0[i] <- gpi_update(alpha[i-1,])
    chi[i,] <- chi_update(gamma = gamma[i-1,], AMAT = AMAT, a = a[i-1,])
    a[i,] <- a_update(AMAT = AMAT, chi = chi[i,])
    pi1[i,] <- spi_update(AMAT = AMAT, a = a[i,])
    pi2[i] <- tpi_update(omega[i-1,])
    alpha[i,] <- alpha_update(X = X, Y = Y, b = b[,,i], group = group, pi0 = pi0[i], alpha = alpha[i-1,], gamma = gamma[i-1,], omega = omega[i-1,], Sigma = Sigma[,,i-1])
    gamma[i,] <- gamma_update(X = X, Y = Y, b = b[,,i], group = group, pi1 = pi1[i,], alpha = alpha[i,], gamma = gamma[i-1,], omega = omega[i-1,], Sigma = Sigma[,,i-1])
    omega[i,] <- omega_update(X = X, Y = Y, b = b[,,i], group = group, pi2 = pi2[i], alpha = alpha[i,], gamma = gamma[i,], omega = omega[i-1,], Sigma = Sigma[,,i-1])
    Z[,,i] <- Z_update(X = X, Y = Y, group = group, alpha = alpha[i,], gamma = gamma[i,], omega = omega[i,])
    beta[,,i] <- beta_update(b = b[,,i], Z = Z[,,i])
    Sigma[,,i] <- Sigma_update(X = X, Y = Y, beta = beta[,,i], b = b[,,i], bsigma = bsigma[i])
  }
  lreturn <- list(group_pi = pi0, snp_pi = pi1, trait_pi = pi2, group_indicator = alpha,
                  snp_indicator = gamma, trait_indicator = omega, b = b, beta = beta,
                  chi = chi, probit_a = a, Sigma = Sigma, biter = c(rep(0, burnin), rep(1, niter)),
                  burnin = burnin, niter = niter,
                  X = X, Y = Y, group = group, bsigma = bsigma, Z = Z)
  return(lreturn)
}
