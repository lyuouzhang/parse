#' @name heatmap_fit
#' @aliases heatmap_fit
#' @title summary plot of globally and pairwise informative variables
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette
#'
#' @description
#' Heatmaps of the data with estimated informative variables and the indicator for pairwise informativeness of each globally informative variables.
#' @usage
#' heatmap_fit(output, y, plot_type = 'info.data', eps.diff = 1e-5,
#' margins = c(5,5), cexRow = 0.5, cexCol = 0.4, lhei = c(0.8,5),
#' lwid=c(0.8,5), adjCol = c(0.8,0.4), sepwidth=c(0.05,0.05))
#' @param output results from parse, apfp, apL1 or nopenalty functions. For the `nopenalty' function, the `short.output' should be FALSE.
#' @param y data.
#' @param plot_type takes two values, 'info.data' or 'info.pair'. 'info.data' is the heatmap of the data with informative variables; 'info.pair' indicates which globally informative variable is pairwise informative for each pair of clusters.
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0. The default value is 1e-5.
#' @param margins parameter in 'heatmap.2' function, 2-dimensional numeric vector containing the margins for column and row names, respectively.
#' @param cexRow parameter in 'heatmap.2' function, positive numbers for the row axis labeling.
#' @param cexCol parameter in 'heatmap.2' function, positive numbers for the column axis labeling.
#' @param lhei parameters in 'heatmap.2' function, visual layout of column height.
#' @param lwid parameters in 'heatmap.2' function, visual layout of column weight.
#' @param adjCol parameters in 'heatmap.2' function, justification of column labels (variables names).
#' @param sepwidth parameters in 'heatmap.2' function, 2-dimensional vector giving the width and height of the separator box
#' @return
#' heatmap of the data with informative variables or heatmap of whether the globally informative variables are pairwise informative for each pair of clusters or not.
#'
#' @references
#' Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2015). gplots: Various R Programming Tools for Plotting Data. R package version 2.17.0.
#' \url{https://CRAN.R-project.org/package=gplots}

#' @examples
#' y <- rbind(matrix(rnorm(120,0,1),ncol=4),
#' matrix(rnorm(120,4,1), ncol=4), matrix(rnorm(120,0,1),ncol=4))
#' output <- parse(K = 3, lambda = 1, y=y)
#' output$mu.hat.best
#' heatmap_fit(output, y, cexRow=1)
#'
#' @keywords external
#' @seealso \link[gplots]{heatmap.2}
#'
#' @export
#'

heatmap_fit <- function(output, y, plot_type='info.data', eps.diff = 1e-5, margins = c(5,5), cexRow = 0.5, cexCol = 0.4, lhei = c(0.8,5),lwid=c(0.8,5), adjCol = c(0.8,0.4), sepwidth=c(0.05,0.05) )
{
  ct.dist = function(x){
    return(sum(dist(x, method ='manhattan') > eps.diff))
  }
  info.index = which(apply(output$mu.hat.best,2,ct.dist) > 0)
  ### only the informative data
  if(plot_type == 'info.data'){
    my_palette <- colorRampPalette(c("white","lightgray", "black"))(n = 256)
    order1 = order(output$s.hat.best)
    heatmap.2(as.matrix(y[order1, info.index]), dendrogram='none',Rowv=F, Colv=T, scale = 'column',density.info = 'none', trace='none', key=F,
              margins = margins, cexRow = cexRow, cexCol = cexCol, col = my_palette,
              lhei = lhei,lwid=lwid, adjCol = adjCol, labRow = order1,
              rowsep = cumsum(table(output$s.hat.best)), sepcolor="white", sepwidth = sepwidth )
  }

  else if (plot_type == 'info.pair'){
    pair_cluster = combn(output$K.best,2)
    pair.info0 = matrix(0, ncol = dim(pair_cluster)[2], nrow = dim(output$mu.hat.best)[2])
    for(j in 1:dim(pair_cluster)[2]){
      temp = which(abs(output$mu.hat.best[pair_cluster[1,j],] - output$mu.hat.best[pair_cluster[2,j],]) > eps.diff)
      pair.info0[temp,j] = 1
    }

    pair.info1 = as.matrix(pair.info0[which(rowSums(pair.info0)!=0),])
    if(dim(pair.info1)[2] > 1){
      colnames(pair.info1) = apply(pair_cluster, 2, function(x){paste('C',x[1],'-C',x[2],sep='')})
      lab.genes = c(1:dim(pair.info1)[1])
      if(length(names(y))>0)
      {
        lab.genes = names(y)[which(rowSums(pair.info0)!=0)]
      }
      black.white <- colorRampPalette(c("white", "black"))(n = 20)
      heatmap.2(as.matrix(pair.info1[order(rowSums(pair.info1)),]), dendrogram='none',Rowv=F, Colv=F, scale = 'none',density.info = 'none', trace='none', key=F,
                margins = margins, cexRow = cexRow, cexCol = cexCol, col = black.white, labRow = lab.genes, srtCol=0, adjCol = adjCol,
                colsep = c(0:dim(pair.info1)[2]),
                rowsep = c(0:dim(pair.info1)[1]),
                sepcolor = "gray",
                sepwidth = sepwidth,
                lhei = lhei,
                lwid = lwid)
    }
  }
}
