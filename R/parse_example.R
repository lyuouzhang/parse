#' @name parse_example
#' @aliases parse_example
#'
#' @title Example for PARSE
#'
#' @import mvtnorm
#' @importFrom mvtnorm rmvnorm
#'
#' @details
#' This is an example of PARSE
#'
#' @examples
#' p = 200
#' mu.all = matrix(0,p,4)
#' mu.all[1:5,1] = 2.5
#' mu.all[1:5,4] = -2.5
#' mu.all[6:10,1] = 1.5
#' mu.all[6:10,2] = 1.5
#' mu.all[6:10,3] = -1.5
#' mu.all[6:10,4] = -1.5
#' Sigma = diag(1,p)
#' set.seed(2023)
#' n = 20
#' y = rbind(rmvnorm(n,mu.all[,1],Sigma),rmvnorm(n,mu.all[,2],Sigma),
#'        rmvnorm(n,mu.all[,3],Sigma),rmvnorm(n,mu.all[,4],Sigma))
#' parse.fit = parse(K = 3:5, lambda = 1:5, y = y, kms.nstart = 10, cores = 8)
#' z.parse = parse.fit$s.hat.best
#'
#' @export

parse_example = function(){
  p = 200
  mu.all = matrix(0,p,4)
  mu.all[1:5,1] = 2.5
  mu.all[1:5,4] = -2.5
  mu.all[6:10,1] = 1.5
  mu.all[6:10,2] = 1.5
  mu.all[6:10,3] = -1.5
  mu.all[6:10,4] = -1.5
  Sigma = diag(1,p)
  n = 20
  set.seed(2023)
  y = rbind(rmvnorm(n,mu.all[,1],Sigma),rmvnorm(n,mu.all[,2],Sigma),rmvnorm(n,mu.all[,3],Sigma),
            rmvnorm(n,mu.all[,4],Sigma))
  parse.fit = parse(K = 3:5, lambda = 2:4, y = y, kms.nstart = 10, cores = 8)
  z.parse = parse.fit$s.hat.best
  mu.parse = parse.fit$mu.hat.best
  return(z.parse=z.parse,mu.parse=mu.parse)
}

