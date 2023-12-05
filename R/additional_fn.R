require(mvtnorm)
## ----- used in sapply to find all the densities
dmvnorm_log <- function(index, mu, sigma, y) {

  ## index = row index of mu
  ## mu = K by p matrix, each row represents one cluster mean
  ## y = n by p data matrix
  ## sigma = p by p covariance matrix (assume same covariance for each cluster)
  ## log.scale = T means output is log of the density

  return(dmvnorm(y, mu[index,], sigma, log=TRUE))
}

## ----- compute the number of unique cluster means for each dimension
## ----- used in compute BIC or GIC
count.mu <- function(mu.j, eps.diff) {
  temp.dist <- as.matrix(dist(mu.j, method = 'manhattan'))
  ct <- length(mu.j[abs(mu.j)>eps.diff])
        ## initial counts (nonzero elements)
  ## --- if exists same means
  temp.dist[upper.tri(temp.dist, diag = T)] <- NA
  if(any(temp.dist < eps.diff, na.rm = T)) {
    temp.index <- which(temp.dist < eps.diff, arr.ind = TRUE)
    temp1 <- mu.j
    ## --- truncated distance so means are not exactly same, make them equal
    for(i in 1:dim(temp.index)[1]){
      temp1[temp.index[i,]] <- min(temp1[temp.index[i,]])
    }
    ct <- length(unique(temp1[abs(temp1)>eps.diff]))
  }
  return(ct)
}



