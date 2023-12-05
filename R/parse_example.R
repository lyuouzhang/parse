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
#' mu.all = matrix(0,30,4)
#' mu.all[1:5,1] = 2.5
#' mu.all[1:5,4] = -2.5
#' mu.all[6:10,1] = 1.5
#' mu.all[6:10,2] = 1.5
#' mu.all[6:10,3] = -1.5
#' mu.all[6:10,4] = -1.5
#' Sigma = diag(4,30)
#' n = 20
#' y = rbind(rmvnorm(n,mu.all[,1],Sigma),rmvnorm(n,mu.all[,2],Sigma),
#'        rmvnorm(n,mu.all[,3],Sigma),rmvnorm(n,mu.all[,4],Sigma))
#' parse.fit = parse(K = 3:5, lambda = 1:5, y = y, kms.nstart = 10, cores = 8)
#' z.parse = parse.fit$s.hat.best
#'
#' @export

parse_example = function(){
  mu.all = matrix(0,30,4)
  mu.all[1:5,1] = 2.5
  mu.all[1:5,4] = -2.5
  mu.all[6:10,1] = 1.5
  mu.all[6:10,2] = 1.5
  mu.all[6:10,3] = -1.5
  mu.all[6:10,4] = -1.5
  Sigma = diag(4,30)
  n = 20
  y = rbind(rmvnorm(n,mu.all[,1],Sigma),rmvnorm(n,mu.all[,2],Sigma),rmvnorm(n,mu.all[,3],Sigma),
            rmvnorm(n,mu.all[,4],Sigma))
  parse.fit = parse(K = 3:5, lambda = 1:5, y = y, kms.nstart = 10, cores = 8)
  z.parse = parse.fit$s.hat.best
  return(z.parse)
}

