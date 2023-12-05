##  ----------------- objective function in stepwise algorithm
##  --- dimension-wise --- input of step.optim
##  --- change pi0 -> alpha (posterior matrix)
obj_parse_fn <- function(mu.j, loc, sigma.j, y.j, alpha, K, lambda, eps.diff = 1e-5){
  ## mu.j = distinct values (in original order)
  ## loc = which clusters have same mean, e.g. c(1,1,2,2,3) means C1=C2, C3=C4
  n = length(y.j)
  if(length(mu.j) == 1){        ## penalty = 0
    mu.j1 <- rep(mu.j, K)       ## re-create mu with length K
    out <- sum(((rep(y.j, times=K) - rep(mu.j1, each=n))^2)*c(alpha))/(2*sigma.j)
  }
  mu.j1 <- mu.j[loc]	    ## same as mu.j%*%loc
  dist.mtx <- dist(mu.j1)
  out <- sum(((rep(y.j, times=K) - rep(mu.j1, each=n))^2)*c(alpha))/(2*sigma.j) + lambda*sum((dist.mtx[dist.mtx > eps.diff])^(-1))    ## truncated differences
  return(out)
}


## ----- stepwise algorithm embeded with optim function with BFGS
parse_backward <- function(j, y, mu.all, alpha, sigma.all, K, lambda, eps.diff = 1e-5, optim.method = 'BFGS'){

  loc1 <- matrix(0,ncol = K, nrow = K)
  obj.value <- c()
  mu.iter <- matrix(0, ncol = K, nrow = K)
  loc1[1,] <- c(1:K)			#starting label

  # input function: obj_parse_fn()
  mu.init1 <- mu.all[,j]
  out0 <- optim(mu.init1[1:K], obj_parse_fn, method = optim.method, loc = loc1[1,], sigma.j = sigma.all[j], y.j = y[,j], alpha = alpha, eps.diff = eps.diff, K = K, lambda = lambda)
  mu.iter[1,]  <- out0$par
  obj.value[1] <- out0$value

  for(q in 1:(K-1)){
    loc.val <- unique(loc1[q,])
    combin1 <- combn(loc.val,2)
    temp.obj <- c()
    temp.mu <- matrix(0, ncol = K, nrow = dim(combin1)[2])
    temp.loc <- matrix(0, ncol = K, nrow = dim(combin1)[2])

    for(i in 1:dim(combin1)[2]){
      temp.loc[i,] <- loc1[q,]
      temp0 <- which(temp.loc[i,] == combin1[2,i])
      temp.loc[i,temp0] <- combin1[1,i]
      temp.loc[i, which(temp.loc[i,] > which(unique(temp.loc[i,]) != c(1:(K-q)))[1])] <- temp.loc[i, which(temp.loc[i,]>which(unique(temp.loc[i,])!=c(1:(K-q)))[1])]-1

      out1 <- optim(mu.init1[1:(K-q)], obj_parse_fn, method = optim.method, loc = temp.loc[i,], sigma.j = sigma.all[j], y.j = y[,j], alpha = alpha, eps.diff = eps.diff, K = K, lambda = lambda)
      temp.obj[i] <- out1$value
      temp.mu[i,] <- out1$par[temp.loc[i,]]
    }
    obj.value[q+1] <- min(temp.obj)
    mu.iter[(q+1),] <- temp.mu[which.min(temp.obj),]
    loc1[(q+1),] <- temp.loc[which.min(temp.obj),]
    if(obj.value[(q+1)] > obj.value[q]){
      break
      }
  }
  return(mu.iter[which.min(obj.value),])
}

