#' @name apL1
#' @aliases apL1
#' @title Model-based Clustering with APL1
#'
#' @importFrom mvtnorm dmvnorm
#' @import stats
#'
#' @description The adaptive \eqn{L_1} penalty was proposed by Pan and Shen (2007). Under the framework of the model-based clustering, APL1 aims to identify the globally informative variables for clustering high-dimensional data.
#'
#' @usage apL1(tuning, K = NULL, lambda = NULL, y, N = 100, kms.iter = 100, kms.nstart = 100,
#'       adapt.kms = FALSE, eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic')
#' apL1(tuning = NULL, K, lambda, y, N = 100, kms.iter = 100, kms.nstart = 100,
#'       adapt.kms = FALSE, eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic')
#'
#' @param tuning A 2-dimensional vector or a matrix with 2 columns, the first column is the number of clusters \eqn{K} and the second column is the tuning parameter \eqn{\lambda} in the penalty term. If this is missing, then \code{K} and \code{lambda} must be provided.
#' @param K The number of clusters \eqn{K}.
#' @param lambda The tuning parameter \eqn{\lambda} in the penalty term.
#' @param y A p-dimensional data matrix. Each row is an observation.
#' @param N The maximum number of iterations in the EM algorithm. The default value is 100.
#' @param kms.iter The maximum number of iterations in kmeans algorithm for generating the starting value for the EM algorithm.
#' @param kms.nstart The number of starting values in K-means.
#' @param adapt.kms A indicator of using the cluster means estimated by K-means to calculate the adaptive parameters in APFP. The default value is FALSE.
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0.
#' @param eps.em The lower bound for the stopping criterion.
#' @param model.crit The criterion used to select the number of clusters \eqn{K}. It is either `bic' for Bayesian Information Criterion or `gic' for Generalized Information Criterion.
#'
#' @details A variable is defined as globally informative if there exists at least one pair of clusters such that \eqn{\mu_{kj} \neq \mu_{k'j}}. Here we assume that each cluster has the same diagonal variance in the model-based clustering. APL1 is in the following form,
#' \deqn{\sum_{j=1}^d \tau_{kj}\sum_{k=1}^K |\mu_{kj}|,}
#' where \eqn{d} is the number of variables in the data, K is the number of clusters, \eqn{\tau_{kj} = \tilde{\mu}_{kj}} is the adaptive parameters. Here we provide two choices for \eqn{\tau_{kj}}. If \code{adapt.kms == TRUE}, \eqn{\tilde{\mu}_{kj}} is the estimates from the K-mean algorithm; otherwise, \eqn{\tilde{\mu}_{kj}} is the estimates from the model-based clustering without penalty.
#'
#' The EM algorithm is used for estimating parameters. Since the EM algorithm depends on the starting values. We use the estimates from K-means with multiple starting points as the starting values.
#'
#' @return This function returns the esimated parameters and some statistics of the optimal model within the given \eqn{K} and \eqn{\lambda}, which is selected by BIC when \code{model.crit = 'bic'} or GIC when \code{model.crit = 'gic'}.
#' \item{mu.hat.best}{The estimated cluster means in the optimal model}
#' \item{sigma.hat.best}{The estimated covariance in the optimal model}
#' \item{p.hat.best}{The estimated cluster proportions in the optimal model}
#' \item{s.hat.best}{The clustering assignments using the optimal model}
#' \item{lambda.best}{The value of \eqn{\lambda} that provide the optimal model}
#' \item{K.best}{The value of \eqn{K} that provide the optimal model}
#' \item{llh.best}{The log-likelihood of the optimal model}
#' \item{gic.best}{The GIC of the optimal model}
#' \item{bic.best}{The BIC of the optimal model}
#' \item{ct.mu.best}{The degrees of freedom in the cluster means of the optimal model}
#'
#' @references Pan, W. and Shen, X. (2007). Penalized model-based clustering with application to variable selection. \emph{The Journal of Machine Learning Research} \bold{8}, 1145--1164.
#'
#' @seealso \code{\link{nopenalty}} \code{\link{apfp}} \code{\link{parse}}
#'
#' @examples
#' y <- rbind(matrix(rnorm(100,0,1),ncol=2), matrix(rnorm(100,4,1), ncol=2))
#' output <- apL1(K = c(1:2), lambda = c(0,0.1), y=y)
#' output$mu.hat.best
#'
#' @keywords external
#'
#' @export


apL1 <- function(tuning=NULL, K=NULL, lambda = NULL,  y, N = 100, kms.iter = 100, kms.nstart = 100, adapt.kms = FALSE, eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic'){
  ## tuning.all: a matrix with 2 columns;
  ##  1st column = K (number of clusters), positive integer;
  ##  2nd column = lambda, nonnegative real number.

  ## check invalid inputs
  y = as.matrix(y)
  if (missing(tuning)){
    if(missing(K) | missing(lambda)){
      stop("Require a matrix of tuning parameters for 'tuning' or vectors for 'K' and 'lambda' ")
    }
    else{
      tuning = as.matrix(expand.grid(K, lambda))
      colnames(tuning)=c()
      }
  }

  if(is.vector(tuning) == TRUE){
    tuning = unlist(tuning)
    if(length(tuning) != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if(dim(y)[1] < tuning[1]){
      stop(" The number of clusters 'K' is greater than the number of observations in 'y' ")
    }
    else if( tuning[1] - round(tuning[1]) > .Machine$double.eps^0.5 | tuning[1] <= 0 | tuning[2] < 0){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }

  else if(is.matrix(tuning) == TRUE){
    if(dim(tuning)[2] != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if(any(dim(y)[1] < tuning[,1])){
      stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
    }
    else if(any(tuning[,1] - round(tuning[,1]) > .Machine$double.eps^0.5 | tuning[,1] <= 0 | tuning[,2] < 0)){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }

  else if(is.matrix(tuning) == FALSE){
    stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
  }

  ### --- end of checking

  ## data dimensions
  n = nrow(y)
  d = ncol(y)

  ### --- choose lambda ---
  s.hat =list()	            # clustering labels
  mu.hat = list()						# cluster means
  sigma.hat = list()				# cluster variance
  p.hat = list() 				    # cluster proportion

  ## ------------ BIC and GIC ----
  llh=c()
  ct.mu=c()
  gic=c()
  bic=c()

  for(j.tune in 1:dim(tuning)[1]){

    K1 = tuning[j.tune,1]
    lambda1 = tuning[j.tune,2]

    if(K1 == 1) {
      mu.hat[[j.tune]] <- colMeans(y)
      sigma.hat[[j.tune]] <- apply(y, 2, var)
      p.hat[[j.tune]] <- 1
      s.hat[[j.tune]] <- rep(1,n)
      llh[j.tune]   <- sum(dmvnorm(x = y, mean = mu.hat[[j.tune]], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
      ct.mu[j.tune] <- sum(abs(mu.hat[[j.tune]]) > eps.diff)
      gic[j.tune]   <- -2*llh[j.tune] + log(log(n))*log(d)*(d + ct.mu[j.tune])
      bic[j.tune]   <- -2*llh[j.tune] + log(n)*(d + ct.mu[j.tune])
    }

    else if(lambda1 == 0)
    {
      temp.out <- nopenalty(K=K1, y=y, N=N, kms.iter=kms.iter, kms.nstart=kms.nstart, eps.diff=eps.diff, eps.em=eps.em, model.crit=model.crit, short.output = FALSE)
      mu.hat[[j.tune]]    <- temp.out$mu.hat.best
      sigma.hat[[j.tune]] <- temp.out$sigma.hat.best
      p.hat[[j.tune]]     <- temp.out$p.hat.best
      s.hat[[j.tune]]     <- temp.out$s.hat.best
      llh[j.tune]         <- temp.out$llh.best
      ct.mu[j.tune]       <- temp.out$ct.mu.best
      gic[j.tune]         <- temp.out$gic.best
      bic[j.tune]         <- temp.out$bic.best
    }

    else{
    sigma.iter <-  matrix(0, ncol = d, nrow = N)	## each row is variances (d dimensions)
    p <- matrix(0, ncol = K1, nrow = N)		## each row is clusters proportions
    mu <- array(0, dim = c(N, K1, d))		## N=iteration, K1=each cluster , d=dimension

  ### --- start from kmeans with multiple starting values	-----
    kms1 <- kmeans(y, centers = K1, nstart= kms.nstart, iter.max = kms.iter)
    kms1.class <- kms1$cluster

    mean0.fn <- function(k){
      if(length(y[kms1.class == k, 1]) == 1) {
        ## if there is only one data in a cluster,
        out <- y[kms1.class == k,]
      }
      else {
        out <- colMeans(y[kms1.class == k,])
      }
      return(out)
    }

    mu[1,,] <- t(sapply(1:K1, mean0.fn))
    sigma.iter[1,] <- apply(y, 2, var)
    p[1, ] <- as.numeric(table(kms1.class))/n

    ## use kmeans' results as the adaptive parameter
    if(adapt.kms == FALSE){
      temp.out <- nopenalty(K=K1, y=y, N=N, kms.iter=kms.iter, kms.nstart=kms.nstart, eps.diff=eps.diff, eps.em=eps.em, model.crit=model.crit, short.output = TRUE)
      mu.no.penal <- temp.out$mu.hat.best
    }
    else {
      mu.no.penal <- mu[1,,]
    }

  ## array of posterior probability computed in E-step
    alpha.temp <- array(0, dim=c(N, n, K1))

    for (t in 1:(N-1)) {
    # E-step: compute all the cluster-wise density
      temp.normal <- sapply(c(1:K1), dmvnorm_log, y = y, mu = mu[t,,], sigma = diag(sigma.iter[t,]))

      alpha.fn = function(k, log.dens = temp.normal, p.temp = p[t,]){
        if(p.temp[k] ==0){out.alpha = rep(0,dim(log.dens)[1])}
        else {
          log.mix1 = sweep(log.dens,2, log(p.temp), '+')
          log.ratio1 = sweep(log.mix1, 1, log.mix1[,k], '-')
          out.alpha = 1/rowSums(exp(log.ratio1))
          out.alpha[which(rowSums(exp(log.ratio1)) == 0)] = 0
        }
        return(out.alpha)
      }

      # E-step: update posterior probabilities alpha
      alpha.temp[(t+1),, ] = sapply(c(1:K1), alpha.fn, log.dens = temp.normal, p.temp = p[t,])

    # M-step  part 1: update cluster proportions p_k
      p[(t+1), ] <- colMeans(alpha.temp[(t+1),,])

    # M-step  part 2: update cluster sigma_j, c() read matrix by columns
      sig.fn = function(index, all.combined, alpha) {  #alpha = alpha[(t+1),,]
        combined = all.combined[,index]
        y.j = combined[1:n]		# n observation's l-th dimension
        mu.j = combined[(n+1):length(combined)]			# K means' l-th dimension
        return(sum((rep(y.j, times= length(mu.j)) - rep(mu.j, each=n))^2*c(alpha))/n)
      }
      all.combined = rbind(y, mu[t,,])   # (n+K1) by d matrix
      sigma.iter[(t+1),] = sapply(c(1:d), sig.fn, all.combined=all.combined, alpha=alpha.temp[(t+1),,])

      apL1.optim <- function(j, y, alpha, sigma.all, mu.no.penal, lambda, K){
        out <- c()
        for(k in 1:K){
          temp1 <- alpha[,k]*y[,j]
          temp2 <- 1-exp(log(lambda)+log(sigma.all[j])-log(abs(sum(temp1)))-log(abs(mu.no.penal[k,j])))
          out[k] <- ifelse(temp2>=0, temp2, 0)*sum(temp1)/sum(alpha[,k])
        }
        return(out)
      }
      mu[(t+1),,] <- sapply(1:d, apL1.optim, y = y, alpha = alpha.temp[(t+1),,], sigma.all = sigma.iter[(t+1),], mu.no.penal = mu.no.penal, lambda = lambda1, K = K1)
    #sapply(combine results in columns = cbind)

      if(sum(abs(mu[(t+1),,] - mu[t,,]))/(sum(abs(mu[t,,])) + 1) + sum(abs(sigma.iter[(t+1),]-sigma.iter[t,]))/sum(abs(sigma.iter[t,])) + sum(abs(p[(t+1),] - p[t,]))/sum(p[t,]) < eps.em){
        s.hat[[j.tune]] <- apply(alpha.temp[(t+1),,], 1, which.max)
        mu.hat[[j.tune]] <- mu[(t+1),,]
        sigma.hat[[j.tune]] <- sigma.iter[(t+1),]
        p.hat[[j.tune]] <- p[(t+1),]
        break
      }
      if(t == N-1) {
        s.hat[[j.tune]] <- apply(alpha.temp[(t+1),,], 1, which.max)
        mu.hat[[j.tune]] <- mu[(t+1),,]
        sigma.hat[[j.tune]] <- sigma.iter[(t+1),]
        p.hat[[j.tune]] <- p[(t+1),]
        warning(paste('not converge when lambda = ', lambda1, 'and K = ', K1, sep=''))
      }
    }
  #ends of each estimates

  ## ------------ BIC and GIC ----
    label.s <- sort(unique(s.hat[[j.tune]]))

    ## with empty clusters
    if(length(label.s) < K1){
      llh[j.tune] <- sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]][label.s]))
      for(k in label.s){
        llh[j.tune] <- llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] == k,], mean = mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log = TRUE))
      }
    }
    ## no empty clusters
    else {
      llh[j.tune] <- sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]]))
      for(k in 1:K1){
        llh[j.tune] <- llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] == k,], mean = mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log = TRUE))
      }
    }

    ct.mu[j.tune] <- sum(apply(mu.hat[[j.tune]], 2, count.mu, eps.diff = eps.diff))
    gic[j.tune] <- -2*llh[j.tune] + log(log(n))*log(d)*(length(label.s) - 1 + d + ct.mu[j.tune])
    bic[j.tune] <- -2*llh[j.tune] + log(n)*(length(label.s) - 1 + d + ct.mu[j.tune])
    }
}  # end of multiple lambda

  if(model.crit == 'bic'){
    index.best = which.min(bic)
  }
  else {
    index.best = which.min(gic)
  }

  gic.best = gic[index.best]
  bic.best = bic[index.best]
  llh.best = llh[index.best]
  ct.mu.best = ct.mu[index.best]

  p.hat.best = p.hat[[index.best]]
  s.hat.best = s.hat[[index.best]]
  mu.hat.best = mu.hat[[index.best]]
  sigma.hat.best = sigma.hat[[index.best]]

  K.best = tuning[index.best,1]
  lambda.best = tuning[index.best,2]

  output = structure(list(K.best = K.best, mu.hat.best=mu.hat.best, sigma.hat.best=sigma.hat.best, p.hat.best=p.hat.best, s.hat.best=s.hat.best, lambda.best=lambda.best, gic.best = gic.best, bic.best=bic.best, llh.best = llh.best, ct.mu.best = ct.mu.best), class = 'parse_fit')

  return(output)
}  # end of function


