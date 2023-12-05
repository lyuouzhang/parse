#' @name nopenalty
#' @aliases nopenalty
#' @title Classical Model-based Clustering
#' @description This function estimates the model-based clustering which is under the framework of finite mixture models.
#'
#' @usage nopenalty(K, y, N = 100, kms.iter = 100, kms.nstart = 100,
#'            eps.diff = 1e-5, eps.em = 1e-5,
#'            model.crit = 'gic', short.output = FALSE)
#'
#' @param K  A vector of the number of clusters
#' @param y  A p-dimensional data matrix. Each row is an observation
#' @param N  The maximum number of iterations in the EM algorithm. The default value is 100.
#' @param kms.iter The maximum number of iterations in the K-means algorithm whose outputs are the starting values for the EM algorithm
#' @param kms.nstart The number of starting values in K-means
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0
#' @param eps.em The lower bound for the stopping criterion.
#' @param model.crit The criterion used to select the number of clusters \eqn{K}. It is either `bic' for Bayesian Information Criterion or `gic' for Generalized Information Criterion.
#' @param short.output A short version of output is needed or not. A short version is used for computing the adaptive parameters in APFP or APL1 methods. The default value is FALSE.
#'
#' @details This function estimates parameters \eqn{\mu}, \eqn{\Sigma}, \eqn{\pi} and the clustering assignments in the model-based clustering using the mixture model,
#' \deqn{y \sim \sum_{k=1}^K \pi_k f(y|\mu_k, \Sigma)}
#' where \eqn{f(y|\mu_k, \Sigma_k)} is the density function of Normal distribution with mean \eqn{\mu_k} and variance \eqn{\Sigma}. Here we assume that each cluster has the same diagonal variance.
#'
#' This function is also used to compute the adaptive parameters for functions \code{\link{apfp}} and \code{\link{apL1}}.
#'
#' @return This function returns the esimated parameters and some statistics of the optimal model within the given \eqn{K} and \eqn{\lambda}, which is selected by BIC when \code{model.crit = 'bic'} or GIC when \code{model.crit = 'gic'}.
#' \item{mu.hat.best}{The estimated cluster means.}
#' \item{sigma.hat.best}{The estimated covariance.}
#' \item{p.hat.best}{The estimated cluster proportions.}
#' \item{s.hat.best}{The clustering assignments.}
#' \item{K.best}{The value of \eqn{K} that provides the optimal model}
#' \item{llh.best}{The log-likelihood of the optimal model}
#' \item{gic.best}{The GIC of the optimal model}
#' \item{bic.best}{The BIC of the optimal model}
#' \item{ct.mu.best}{The degrees of freedom in the cluster means of the optimal model}
#'
#' @references Fraley, C., & Raftery, A. E. (2002) Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American statistical Association} \bold{97(458)}, 611--631.
#'
#' @seealso \code{\link{apfp}} \code{\link{apL1}} \code{\link{parse}}
#'
#' @examples
#' y <- rbind(matrix(rnorm(100,0,1),ncol=2), matrix(rnorm(100,4,1), ncol=2))
#' output <- nopenalty(K = c(1:2), y)
#' output$mu.hat.best
#'
#' @keywords external
#'
#' @export


nopenalty <- function(K, y, N = 100, kms.iter = 100, kms.nstart = 100, eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic', short.output = FALSE){

  if(missing(K) | any(dim(y)[1] < K)){
    stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
    }

  ## data dimensions
  n = nrow(y)
  D = ncol(y)

  ## outputs
  s.hat =list()	            # clustering labels
  mu.hat = list()						# cluster means
  sigma.hat = list()				# cluster variance
  p.hat = list() 				    # cluster proportion

  ## ------------ BIC and GIC ----
  llh=c()
  ct.mu=c()
  gic=c()
  bic=c()

  for(j.tune in 1:length(K)){

    if(K[j.tune] == 1) {
      mu.hat[[j.tune]] <- colMeans(y)
      sigma.hat[[j.tune]] <- apply(y, 2, var)
      p.hat[[j.tune]] <- 1
      s.hat[[j.tune]] <- rep(1,n)
      llh[j.tune] <- sum(dmvnorm(x = y, mean = mu.hat[[j.tune]], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
      ct.mu[j.tune] <- sum(abs(mu.hat[[j.tune]]) > eps.diff)
      gic[j.tune] <- -2*llh[j.tune] + log(log(n))*log(D)*(D + ct.mu[j.tune])
      bic[j.tune] <- -2*llh[j.tune] + log(n)*(D + ct.mu[j.tune])
    }
    else{
      sigma.iter =  matrix(0, ncol = D, nrow = N)		#each row is variances (D dimensions)
      p = matrix(0, ncol = K[j.tune], nrow = N)		# each row is clusters proportions
      mu = array(0, dim = c(N, K[j.tune], D))		# N=iteration, K=each cluster , D=dimension

  ### --- start from kmeans with multiple starting values	-----
      kms1 = kmeans(y, centers = K[j.tune], nstart = kms.nstart, iter.max = kms.iter)
      kms1.class = kms1$cluster

      mean0.fn = function(k){
        if(length(y[kms1.class == k,1]) == 1) {out = y[kms1.class == k,]}
        else {	out = colMeans(y[kms1.class == k,])}
        return(out)
      }

      mu[1,,] = t(sapply(1:K[j.tune], mean0.fn))
      sigma.iter[1,] = apply(y, 2, var)
      p[1, ] = as.numeric(table(kms1.class))/n

      ## array of posterior probability computed in E-step
      alpha.temp = array(0, dim = c(N, n, K[j.tune]))

      for (t in 1:(N-1)) {
        # E-step: compute all the cluster-wise density
        temp.normal = sapply(c(1:K[j.tune]), dmvnorm_log, y = y, mu = mu[t,,], sigma = diag(sigma.iter[t,]))

        alpha.fn = function(k, log.dens = temp.normal, p.temp = p[t,]){
          if(p.temp[k] == 0){
            out.alpha = rep(0,dim(log.dens)[1])
            }
          else {
            log.mix1 = sweep(log.dens,2, log(p.temp), '+')
            log.ratio1 = sweep(log.mix1, 1, log.mix1[,k], '-')
            out.alpha = 1/rowSums(exp(log.ratio1))
            out.alpha[which(rowSums(exp(log.ratio1)) == 0)] = 0
          }
          return(out.alpha)
        }

        # E-step: update posterior probabilities alpha
        alpha.temp[(t+1),, ] = sapply(c(1:K[j.tune]), alpha.fn, log.dens = temp.normal, p.temp = p[t,])

        # M-step  part 1: update cluster proportions p_k
        p[(t+1), ] = colMeans(alpha.temp[(t+1),,])

        # M-step  part 2: update cluster sigma_j, c() read matrix by columns
        sig.est.fn = function(index, all.combined, alpha) {
          combined = all.combined[,index]
          y.l = combined[1:n]		# n observation's l-th dimension
          mu.l = combined[(n+1):length(combined)]			# K[j.tune] means' l-th dimension
          return(sum((rep(y.l, times=K[j.tune]) - rep(mu.l, each=n))^2*c(alpha))/n)
        }
        all.combined = rbind(y, mu[t,,])   # (n+K[j.tune]) by D matrix
        sigma.iter[(t+1),] = sapply(c(1:D), sig.est.fn, all.combined=all.combined, alpha = alpha.temp[(t+1),,])

        # M-step  part 3: update cluster means mu_j
        mu.est.fn = function(k.index, y, alpha){
          mu.k1 = alpha[,k.index] %*% y/sum(alpha[,k.index])
          return(mu.k1)
        }
        mu[(t+1),,] = t(sapply(1:K[j.tune], mu.est.fn, y = y, alpha = alpha.temp[(t+1),,]))
        #sapply(combine results in columns = cbind)

        # use relative difference as stopping criterion
        if(sum(abs(mu[(t+1),,] - mu[t,,]))/(sum(abs(mu[t,,])) + 1)+ sum(abs(sigma.iter[(t+1),]-sigma.iter[t,]))/sum(abs(sigma.iter[t,])) + sum(abs(p[(t+1),] - p[t,]))/sum(p[t,]) < eps.em) {
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,], 1, which.max)
          mu.hat[[j.tune]] = mu[(t+1),,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[(t+1),]					# cluster variance
          p.hat[[j.tune]] = p[(t+1),]
          break
        }

        if(t==N-1) {
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,], 1, which.max)
          mu.hat[[j.tune]] = mu[(t+1),,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[(t+1),]					# cluster variance
          p.hat[[j.tune]] = p[(t+1),]
          warning(paste('not converge when lambda = 0 and K=',K[j.tune],sep=''))
        }
      }

      ## ------------ BIC and GIC ----
      label.s = sort(unique(s.hat[[j.tune]]))

      ## with empty clusters
      if(length(label.s) < K[j.tune]){
        llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]][label.s]))
        for(k in label.s){
          llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
        }
      }
      ## no empty clusters
      else {
        llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]]))
        for(k in 1:K[j.tune]){
          llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
        }
      }
      ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, eps.diff = eps.diff))
      gic[j.tune] = -2*llh[j.tune] + log(log(n))*log(D)*(length(label.s)-1+D+ct.mu[j.tune])
      bic[j.tune] = -2*llh[j.tune] + log(n)*(length(label.s)-1+D+ct.mu[j.tune])
    }
  }
 # end of multiple K

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

  K.best = K[index.best]

  if(short.output == TRUE ){
    output = list(mu.hat.best = mu.hat.best)
  }
  else {
    output = structure(list(K.best = K.best, mu.hat.best=mu.hat.best, sigma.hat.best=sigma.hat.best, p.hat.best=p.hat.best, s.hat.best=s.hat.best, gic.best = gic.best, bic.best=bic.best, llh.best = llh.best, ct.mu.best = ct.mu.best), class = 'parse_fit')
  }
  return(output)
}
