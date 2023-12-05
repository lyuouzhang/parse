#' @name summary
#' @aliases summary.parse_fit
#' @title summary of the clustering results
#'
#' @description
#' Summary of the globally informative variables
#' @usage
#' summary(output, y, eps.diff = 1e-5)
#'
#' @param output results from parse, apfp, apL1 or nopenalty functions. For the `nopenalty' function, the `short.output' should be FALSE.
#' @param y data.
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0.
#'
#' @return
#' \item{num.info}{the number of globally informative variables}
#' \item{perc.info}{the percentage of globally informative variables}
#' \item{info.name}{the variable names of the globally informative variables, if the data have no variable name, 'info.name' is the index of the variable}
#' @examples
#' y <- rbind(matrix(rnorm(120,0,1),ncol=4),
#' matrix(rnorm(120,4,1), ncol=4), matrix(rnorm(120,0,1),ncol=4))
#' output <- parse(K = c(1:2), lambda = c(0,1), y=y)
#' output$mu.hat.best
#' summary(output, y)
#'
#' @export
#'

summary <- function(output, y, eps.diff = 1e-5)
{
  UseMethod("summary",output)
}
#' @export
summary.parse_fit <- function(output, y, eps.diff = 1e-5)
{
  ct.dist = function(x){
    return(sum(dist(x, method ='manhattan') > eps.diff))
  }

  pair.diff = function(mu){
    n.diff = sum(apply(mu,2,ct.dist)>0)
    return(n.diff = n.diff)
  }

  num.info = pair.diff(output$mu.hat.best)
  perc.info = pair.diff(output$mu.hat.best)/dim(output$mu.hat.best)[2]
  info.index = which(apply(output$mu.hat.best,2,ct.dist) > 0)
  if(length(names(y)) > 0)
  {
    info.name = names(y)[info.index]
    output = list(num.info = num.info, perc.info = perc.info, info.name = info.name)
  }

  else {
    output = list(num.info = num.info, perc.info = perc.info, info.name = info.index)
  }
  return(output)
}
