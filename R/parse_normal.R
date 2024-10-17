#' @name parse_normal
#' @aliases parse_normal
#' @title Model-based Clustering with PARSE
#'
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @import utils
#' @description The PAirwise Reciprocal fuSE (PARSE) penalty was proposed by Wang, Zhou and Hoeting (2016). Under the framework of the model-based clustering, PARSE aims to identify the pairwise informative variables for clustering, especially for high-dimensional data.
#'
#' @usage parse(tuning, K = NULL, lambda = NULL, y, N = 100, kms.iter = 100, kms.nstart = 100,
#'       eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic', backward = TRUE, cores=2)
#' parse(tuning = NULL, K, lambda, y, N = 100, kms.iter = 100, kms.nstart = 100,
#'       eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic', backward = TRUE, cores=2)
#'
#' @param tuning A 2-dimensional vector or a matrix with 2 columns, the first column is the number of clusters \eqn{K} and the second column is the tuning parameter \eqn{\lambda} in the penalty term. If this is missing, then \code{K} and \code{lambda} must be provided.
#' @param K The number of clusters \eqn{K}.
#' @param lambda The tuning parameter \eqn{\lambda} in the penalty term.
#' @param y A p-dimensional data matrix. Each row is an observation.
#' @param N The maximum number of iterations in the EM algorithm. The default value is 100.
#' @param kms.iter The maximum number of iterations in kmeans algorithm for generating the starting value for the EM algorithm.
#' @param kms.nstart The number of starting values in K-means.
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0.
#' @param eps.em The lower bound for the stopping criterion.
#' @param model.crit The criterion used to select the number of clusters \eqn{K}. It is either `bic' for Bayesian Information Criterion or `gic' for Generalized Information Criterion.
#' @param backward Use the backward selection algorithm when it equals to "TRUE", otherwise select all the possible subsets.
#' @param cores The number of cores which can be used in parallel computing.
#' @param min.iter The number of iteration we start to remove non-informative features.
#' @details The j-th variable is defined as pairwise informative for a pair of clusters \eqn{C_k} and \eqn{C_{k'}} if \eqn{\mu_{kj} \neq \mu_{k'j}}. Also, a variable is globally informative if it is pairwise informative for at least one pair of clusters. Here we assume that each cluster has the same diagonal variance in the model-based clustering. PARSE is in the following form,
#' \deqn{\sum_{j=1}^{d}\sum_{k<k'}|\mu_{kj} - \mu_{k'j}|^{-1} \mathbf{I}(\mu_{kj} \neq \mu_{k'j}).}
#' where \eqn{d} is the number of variables in the data.
#'
#' The estimation uses the backward searching algorithm embedded in the EM algorithm. Since the EM algorithm depends on the starting values. We use the estimates from K-means with multiple starting points as the starting values. Please check the paper for details of the algorithm. In this function we use parallel computing to estimate cluster means for each dimension. The default number of cores to be used is 2, which can be specified by users.
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
#' @references Wang, L., Zhou, W. and Hoeting, J. (2016) Identification of Pairwise Informative Features for Clustering Data. \emph{preprint}.
#'
#' @seealso \code{\link[stats]{optim}} \code{\link{nopenalty}} \code{\link{apL1}} \code{\link{apfp}} \code{\link[foreach]{foreach}} \code{\link{doParallel}}
#'
#' @examples
#' y <- rbind(matrix(rnorm(120,0,1),ncol=3), matrix(rnorm(120,4,1), ncol=3))
#' output <- parse(K = c(1:2), lambda = c(0,1), y=y, cores=2)
#' output$mu.hat.best
#' @keywords external
#'
#' @export

utils::globalVariables("j")
parse <- function(tuning=NULL, K=NULL, lambda = NULL, y, N = 100, kms.iter = 100, kms.nstart = 10,
                  eps.diff = 1e-5, eps.em = 1e-5, model.crit = 'gic', backward = TRUE, cores = 2, min.iter = 1,
                  pca_adjust = 0.1){
  ## tuning: a matrix with 2 columns;
  ##  1st column = K (number of clusters), positive integer;
  ##  2nd column = lambda, nonnegative real number.

  ## ----- check invalid numbers
  y = as.matrix(y)
  if(missing(tuning)){
    if(missing(K)||missing(lambda)){
      stop("Require a matrix of tuning parameters for 'tuning' or vectors for 'K' and 'lambda' ")
    }else{
      tuning = as.matrix(expand.grid(K, lambda))
      colnames(tuning)=c()
    }
  }

  if(is.vector(tuning) == TRUE){
    tuning = unlist(tuning)
    if(length(tuning) != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }else if(dim(y)[1] < tuning[1]){
      stop(" The number of clusters 'K' is greater than the number of observations in 'y' ")
    }else if( tuning[1] - round(tuning[1]) > .Machine$double.eps^0.5 | tuning[1] <= 0 | tuning[2] < 0){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }else if(is.matrix(tuning) == TRUE){
    if(dim(tuning)[2] != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }else if(any(dim(y)[1] < tuning[,1])){
      stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
    }else if(any(tuning[,1] - round(tuning[,1]) > .Machine$double.eps^0.5 | tuning[,1] <= 0 | tuning[,2] < 0)){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }else if(is.matrix(tuning) == FALSE){
    stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
  }

  if(detectCores() < cores){
    stop("The number of cores is too large.")
    }else{
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  }

  ## data dimensions
  n = nrow(y)
  d = ncol(y)

  ## PCA adjustment
  ## default 10% informative features
  if(!is.null(pca_adjust)){
    svd.result = svd(y)
    parse.positive = sort(order(apply(abs(svd.result$v[,1:4]),1,sum),decreasing = TRUE)[1:floor(d*pca_adjust)])
  }else{
    parse.positive = 1:d
  }


  ### outputs
  s.hat =list()	            # clustering labels
  mu.hat = list()						# cluster means
  sigma.hat = list()				# cluster variance
  p.hat = list() 				    # cluster proportion
  feature.hat = list()      # working features

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
    }else if(lambda1 == 0){
      temp.out <- nopenalty(K=K1, y=y, N=N, kms.iter=kms.iter, kms.nstart=kms.nstart, eps.diff=eps.diff, eps.em=eps.em, model.crit=model.crit, short.output = FALSE)
      mu.hat[[j.tune]]    <- temp.out$mu.hat.best
      sigma.hat[[j.tune]] <- temp.out$sigma.hat.best
      p.hat[[j.tune]]     <- temp.out$p.hat.best
      s.hat[[j.tune]]     <- temp.out$s.hat.best
      llh[j.tune]         <- temp.out$llh.best
      ct.mu[j.tune]       <- temp.out$ct.mu.best
      gic[j.tune]         <- temp.out$gic.best
      bic[j.tune]         <- temp.out$bic.best
    }else{
      #sigma.iter <-  array(0, dim = c(N,d,d))	## first index is for iteration, the other two for covariance (d by d)
      sigma.iter = list()
      for(i in 3:4){
        #sigma.iter[[i]] = bigmemory::big.matrix(d,d)
        sigma.iter[[i]] = 1:d
      }
      p <- matrix(0, ncol = K1, nrow = N)		## each row is clusters proportions
      mu <- array(0, dim = c(N, K1, d))		## N=iteration, K1=each cluster , d=dimension

      ### --- start from spectral clustering-----
      #if(d<=2000){
      #  kms1 = kmeans(y, centers = K1, nstart = kms.nstart, iter.max = kms.iter)
      #}else{
      #  y0 = y[,order(apply(y,2,var),decreasing = TRUE)[1:2000]]
      #  kms1 = kmeans(y0, centers = K1, nstart = kms.nstart, iter.max = kms.iter)
      #}
      kms1 = kernlab::specc(y[,parse.positive], centers = K1)
      kms1.class = as.numeric(kms1)

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

      #sigmahat with sample covariance matrix
      sigma.est.fn <- function(y, mu, alpha ,n, d, K){
        sigma = matrix(0,d,d)
        for(k in seq(K)){
          v.tmp = t( t(y) - mu[k,] )
          sigma <- sigma + t(v.tmp*alpha[,k])%*%v.tmp/n
          # for(i in seq(n)){
          #   v.tmp <- y[i,] -mu[k,]
          #   sigma <- sigma + alpha[i,k] * tcrossprod(v.tmp)/n
          # }
        }
        return(sigma)
      }

      #sigmahat with for only diagonal
      sigma.est.fn.single <- function(y, mu, alpha ,n, d, K){
        sigma = rep(0,d)
        for(k in seq(K)){
          v.tmp = t( t(y) - mu[k,] )
          sigma = sigma + apply(t(v.tmp^2)*alpha[,k],1,sum)/n
          # for(i in seq(n)){
          #   for(j in seq(d)){
          #     v.tmp <- y[i,j] -mu[k,j]
          #     sigma[j] <- sigma[j] + alpha[i,k] * v.tmp^2/n
          #   }
          # }
        }
        return(sigma)
      }

      mu[1,,] <- t(sapply(1:K1, mean0.fn))
      p[1, ] <- as.numeric(table(kms1.class))/n
      #sigma.iter[1,,] <- sigma.est.fn(y,mu[1,,],matrix(p[1, ],n,K1,byrow = TRUE),n, d, K1)
      sigma.iter[[1]] = sigma.iter[[2]] = apply(y,2,var)

      ## array of posterior probability computed in E-step
      alpha.temp <- array(0, dim=c(N,n,K1))

      ### difference between means
      mu.diff.create = function(K){
        A = c()
        for(k in 1:K){
          A = rbind(A, cbind(matrix(0,K-k,k-1),rep(1,K-k),-diag(K-k)))
        }
        return(A)
      }

      alpha.fn.multi <- function(y.tmp, mu.tmp, sigma.tmp, pro.tmp, n, K, seq_K = 2:K){
        ret <- matrix(0,n,K)
        gamma.tmp = t(solve(sigma.tmp,mu.tmp[,1]-mu.tmp[,-1]))
        for(k in 2:K){
          ret[,k] <- -gamma.tmp[k-1,]%*%(y.tmp-(mu.tmp[,1]+mu.tmp[,k])/2)
        }
        exp.ret= exp(ret + t(matrix(log(pro.tmp),K,n)))
        exp.ret[which(exp.ret==Inf)] = 1e6
        ret <- exp.ret/apply(exp.ret,1,sum)
        return(ret)
      }

      alpha.fn.single <- function(y.tmp, mu.tmp, sigma.tmp, pro.tmp, n, K, seq_K = 2:K){
        ret <- matrix(0,n,K)
        gamma.tmp = t((mu.tmp[,1]-mu.tmp[,-1])/sigma.tmp)
        for(k in 2:K){
          ret[,k] <- -gamma.tmp[k-1,]%*%(y.tmp-(mu.tmp[,1]+mu.tmp[,k])/2)
        }
        exp.ret= exp(ret + t(matrix(log(pro.tmp),K,n)))
        exp.ret[which(exp.ret==Inf)] = 1e6
        ret <- exp.ret/apply(exp.ret,1,sum)
        return(ret)
      }

      mu_diff = mu.diff.create(K1)

      for (t in 1:(N-1)) {
        # E-step: compute all the cluster-wise density
        # only involve informative variables
        # start from min.iter
        if(t>min.iter){
          parse.diff = mu_diff%*%mu[t,,]
          parse.positive = which(apply(abs(parse.diff),2,sum)>1e-4)
        }

        parse.positive = intersect(parse.positive,sigma.iter[[3]])
        parse.position = which(sigma.iter[[3]] %in% parse.positive)

        if(length(parse.positive)<=2){
          s.hat[[j.tune]] = apply(alpha.temp[t,,],1,which.min)		# clustering labels
          mu.hat[[j.tune]] = mu[t,,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[[2]]				# cluster variance
          p.hat[[j.tune]] = p[t,]
          break
        }


        # E-step: update posterior probabilities pi0
        if(length(parse.positive)<n-2){
          if(!is.matrix(sigma.iter[[1]])){
            sigma.tmp = diag(sigma.iter[[1]][parse.position])
          }else{
            sigma.tmp = sigma.iter[[1]][parse.position,parse.position]
          }
          alpha.temp[(t+1),, ] <- alpha.fn.multi(y.tmp=t(y[,parse.positive]), mu.tmp=t(mu[t,,parse.positive]),
                                           sigma.tmp = sigma.tmp,pro.tmp = p[t,], n, K1)
        }else{
          if(is.matrix(sigma.iter[[1]])){
            sigma.tmp = diag(sigma.iter[[1]][parse.position,parse.position])
          }else{
            sigma.tmp = sigma.iter[[1]][parse.position]
          }
          alpha.temp[(t+1),, ] <- alpha.fn.single(y.tmp=t(y[,parse.positive]), mu.tmp=t(mu[t,,parse.positive]),
                                                 sigma.tmp = sigma.tmp,pro.tmp = p[t,], n, K1)
        }

        # M-step  part 1: update cluster proportions p_k

        p[(t+1), ] <- colMeans(alpha.temp[(t+1),,])

        # M-step  part 2: update cluster sigma_j, c() read matrix by columns
        if(length(parse.positive)<n-2){
          sigma.iter[[2]] <- sigma.est.fn(y[,parse.positive],mu[t,,parse.positive],alpha.temp[(t+1),,],n, length(parse.positive), K1)
          sigma.iter[[4]] = parse.positive
        }else{
          sigma.iter[[2]] <- sigma.est.fn.single(y[,parse.positive],mu[t,,parse.positive],alpha.temp[(t+1),,],n, length(parse.positive), K1)
          sigma.iter[[4]] = parse.positive
        }

        #mu[(t+1),,] <- foreach(j=1:d, .combine=cbind) %dopar% {parse_backward(j=j, y=y, mu.all = mu[t,,], alpha=alpha.temp[(t+1),,], sigma.all = diag(sigma.iter[[2]]), eps.diff = eps.diff, K=K1, lambda = lambda1, optim.method ='BFGS')}
        #mu[(t+1),,] <- foreach(j=1:d, .combine=cbind) %dopar% {parse_backward_single(y=y[,j], mu.all = mu[t,,j], alpha=alpha.temp[(t+1),,], eps.diff = eps.diff, K=K1, lambda = lambda1, optim.method ='BFGS')}

        # update informative features and 10% non-informative features
        mu[(t+1),,] = mu[t,,]
        parse.update = sort(c(parse.positive,sample(setdiff(1:d,parse.positive),floor((d-length(parse.positive))*0.1))))
        mu[(t+1),,parse.update] <- foreach(j=parse.update, .combine=cbind) %dopar% {parse_backward_single_neighbor(y=y[,j], mu.all = mu[t,,j], alpha=alpha.temp[(t+1),,], eps.diff = eps.diff, K=K1, lambda = lambda1, optim.method ='BFGS')}


        if(sum(abs(mu[(t+1),,]-mu[t,,]))/(sum(abs(mu[t,,]))+1) + sum(abs(p[(t+1),] - p[t,]))/sum(p[t,]) < eps.em) {
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,], 1, which.min)
          mu.hat[[j.tune]] = mu[(t+1),,]				# cluster means
          sigma.hat[[j.tune]] = sigma.iter[[2]]		# cluster variance
          p.hat[[j.tune]] = p[(t+1),]
          feature.hat[[j.tune]] = sigma.iter[[4]]
          break
        }else{
          sigma.iter[[1]] = sigma.iter[[2]]
          sigma.iter[[3]] = sigma.iter[[4]]
        }

        if(t==N-1) {
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,],1,which.min)		# clustering labels
          mu.hat[[j.tune]] = mu[(t+1),,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[[2]]				# cluster variance
          p.hat[[j.tune]] = p[(t+1),]
          feature.hat[[j.tune]] = sigma.iter[[4]]
          warning(paste('not converge when lambda = ', lambda1, ' and K = ', K1, sep=''))
        }
      }

      ## ------------ BIC and GIC ----
      label.s = sort(unique(s.hat[[j.tune]]))

      if(length(feature.hat[[j.tune]])<n-2){
        if(!is.matrix(sigma.hat[[j.tune]])){
          sigma.tmp = diag(sigma.hat[[j.tune]])
        }else{
          sigma.tmp = sigma.hat[[j.tune]]
        }
        ## with empty clusters
        if(length(label.s) < K1){
          llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]][label.s]))
          for(k in label.s){
            llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = sigma.tmp, log=TRUE))
          }
        }else {
          ## no empty clusters
          llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]]))
          for(k in 1:K1){
            llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = sigma.tmp, log=TRUE))
          }
        }
        ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, eps.diff=eps.diff))
        gic[j.tune] = -2*llh[j.tune] + log(log(n))*log(d)*(length(label.s)-1+d+ct.mu[j.tune])
        bic[j.tune] = -2*llh[j.tune] + log(n)*(length(label.s)-1+d+ct.mu[j.tune])
      }else{
        if(is.matrix(sigma.hat[[j.tune]])){
          sigma.tmp = diag(sigma.hat[[j.tune]])
        }else{
          sigma.tmp = sigma.hat[[j.tune]]
        }
        ## with empty clusters
        if(length(label.s) < K1){
          llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]][label.s]))
          for(k in label.s){
            llh[j.tune] = llh[j.tune] + sum(dnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sd = sqrt(sigma.tmp), log=TRUE))
          }
        }else {
          ## no empty clusters
          llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]]))
          for(k in 1:K1){
            llh[j.tune] = llh[j.tune] + sum(dnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sd = sqrt(sigma.tmp), log=TRUE))
          }
        }
        ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, eps.diff=eps.diff))
        gic[j.tune] = -2*llh[j.tune] + log(log(n))*log(d)*(length(label.s)-1+d+ct.mu[j.tune])
        bic[j.tune] = -2*llh[j.tune] + log(n)*(length(label.s)-1+d+ct.mu[j.tune])
      }
    }
    #ends of one K, lambda
  }
  # end of multiple lambda

  if(model.crit == 'bic'){
    index.best = which.min(bic)
  }else {
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
  feature.hat.best = feature.hat[[index.best]]

  K.best = tuning[index.best,1]
  lambda.best = tuning[index.best,2]

  output = structure(list(K.best = K.best, mu.hat.best=mu.hat.best, sigma.hat.best=sigma.hat.best, p.hat.best=p.hat.best, s.hat.best=s.hat.best, lambda.best=lambda.best, feature.hat.best=feature.hat.best, gic.best = gic.best, bic.best=bic.best, llh.best = llh.best, ct.mu.best = ct.mu.best), class = 'parse_fit')
  stopCluster(cl)
  return(output)
}

