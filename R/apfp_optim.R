#' @keywords internal
## ------------ dimension-wise, based on local quadratic qpproximation
apfp_LQA <- function(j, y, mu.t.all, mu.no.penal, sigma.all, alpha, lambda, iter.num = 20, eps.LQA = 1e-5, eps.diff = 1e-5){

  ## mu.t.all : K by p mean matrix from previous EM-step
  ## mu.no.penal : K by p mean matrix of unpenalized estimates (lambda=0)
  ## y: n by p data matrix
  ## sigma.all: p by p diagnal covariance matrix
  ## alpha: n by K posterior probability matrix
  ## iter.num: max iterations in local quadratic approximation (LQA)
  ## eps.LQA: LQA stop criterion
  ## eps.diff: lower bound of mean difference
  ## lambda: tuning parameter

  mu.t.j <- mu.t.all[,j]
  mu.no.j <- mu.no.penal[,j]
  y.j <- y[,j]
  sigma.j <- sigma.all[j]

  n <- length(y.j)
  K <- length(mu.t.j)

  y.j.tilde <- rep(y[,j], times = K)
  A <- diag(c(alpha)/(2*sigma.j))
  X <- diag(K)[rep(1:K, each = n),]

  mat1 <- t(X)%*%A%*%X
  vec1 <- t(X)%*%A%*%y.j.tilde

  ## if distance < eps.diff, make them equal.
  dist.mu.no.j <- dist(mu.no.j, method ='manhattan')
  dist.mu.no.j[dist.mu.no.j < eps.diff] <- eps.diff
  log.tau <- c(-log(dist(mu.no.j, method ='manhattan')))

  mu.j.hat <- matrix(0, nrow=iter.num+1,ncol=K)
  mu.j.hat[1,] <- mu.t.j

  for(s in 1:iter.num){
    dist.old.tmp0 <- as.matrix(dist(mu.j.hat[s,], method = 'manhattan'))
    dist.old.tmp0[upper.tri(dist.old.tmp0, diag = TRUE)] <- NA
    index0 <- which(is.na(dist.old.tmp0) == FALSE, arr.ind = TRUE)
    dist.old.tmp1 <- dist.old.tmp0[lower.tri(dist.old.tmp0, diag = FALSE)]
    dist.old.tmp1[dist.old.tmp1 < eps.diff] <- eps.diff

    c.kk <- sqrt(exp(log.tau - log(2*dist.old.tmp1)))
    tmp.fn1 <- function(s.i){
      vec.tmp1 <- rep(0,K)
      vec.tmp1[index0[s.i,]] <- c(-1,1)
      return(vec.tmp1)
    }

    if(dim(index0)[1]==1){
      combn.mat <- tmp.fn1(1)
      B <- combn.mat%*%t(combn.mat)
    }
    else {
      combn.mat <- t( sapply(1:dim(index0)[1], tmp.fn1) )
      B <- t(combn.mat)%*%t(diag(c.kk))%*%diag(c.kk)%*%combn.mat
    }
    mu.j.hat[(s+1),] <- solve(mat1+lambda*B)%*%vec1

    ## stopping criterion = relative distance (+1) needed when the vector is 0.
    ## for computing stability, let small value to be exact zero
    if(sum(abs(mu.j.hat[(s+1),] - mu.j.hat[s,]))/(sum(abs(mu.j.hat[s,])) + 1) < eps.LQA){
      out <- mu.j.hat[(s+1),]
      out[which(abs(out) < eps.diff)] <- 0
      break
    }

    else if(s == iter.num){
      out <- mu.j.hat[(s+1),]
      out[which(abs(out) < eps.diff)] <- 0
      warning(paste('LQA for estimating mu does not converge for K=', K, 'lambda=', lambda, sep=''))
    }
  }
  return(out)
}

