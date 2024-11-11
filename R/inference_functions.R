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

bfdr <- function(pips, cutoff, tol = 0.01) {
  phat <- 1 - pips
  tseq <- seq(min(phat), 1, by=1e-4)
  bfdrt <- numeric(length(tseq))
  counter <- 1
  for(t in tseq) {
    dk <- as.numeric(phat < t)
    bfdrt[counter] <- sum(phat * dk) / sum(dk)
    counter <- counter + 1
  }
  diffs <- abs(cutoff - bfdrt)
  tchoice <- min(which(diffs == min(diffs, na.rm = TRUE)))
  fdr_est <- bfdrt[tchoice]
  to_return <- c(1-tseq[tchoice], fdr_est)
  names(to_return) <- c("cutoff value", "estimated bfdr")
  #return(1 - tseq[tchoice])
  if(diffs[tchoice]>tol) {
    warning("tolerance exceeded: check results")
  }
  return(to_return)
}
