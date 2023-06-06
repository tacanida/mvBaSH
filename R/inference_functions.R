#inference functions

get_post_med <- function(mcmc_out) {
  beta <- mcmc_out$beta
  burnin <- mcmc_out$burnin
  return(apply(beta[,,-c(1:burnin)], c(1,2), median))
}

get_overall_pip <- function(mcmc_out) {
  beta <- mcmc_out$beta
  burnin <- mcmc_out$burnin
  beta <- beta[,,-c(1:burnin)]
  return(rowMeans(apply(beta, 3, FUN = function(x) rowSums(x) != 0)))
}

