#inference functions

#get posterior median
pos_med <- function(output) {
  return(apply(output$beta[,,-c(1:output$burnin)], c(1,2), median))
}

#get overall pip
get_pips <- function(output) {
  indsums <- apply(output$Z[,,-c(1:output$burnin)], 3, FUN = function(x) rowSums(x))
  return(apply(indsums, 1, FUN = function(x) mean(x>0)))
}

#find fdr cutoff
bfdr <- function(pips, cutoff) {
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
  return(1 - tseq[tchoice])
}

#this one doesn't really work...
bfdr2 <- function(pips, cutoff) {
  phat <- 1 - pips
  tseq <- unique(phat)
  bfdrt <- numeric(length(tseq))
  counter <- 1
  for(t in tseq) {
    dk <- as.numeric(phat < t)
    bfdrt[counter] <- sum(phat * dk) / sum(dk)
    counter <- counter + 1
  }
  diffs <- abs(cutoff - bfdrt)
  tchoice <- min(which(diffs == min(diffs, na.rm = TRUE)))
  return(1 - tseq[tchoice])
}
