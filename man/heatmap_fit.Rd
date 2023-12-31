% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_parse_fit.R
\name{heatmap_fit}
\alias{heatmap_fit}
\title{summary plot of globally and pairwise informative variables}
\usage{
heatmap_fit(output, y, plot_type = 'info.data', eps.diff = 1e-5,
margins = c(5,5), cexRow = 0.5, cexCol = 0.4, lhei = c(0.8,5),
lwid=c(0.8,5), adjCol = c(0.8,0.4), sepwidth=c(0.05,0.05))
}
\arguments{
\item{output}{results from parse, apfp, apL1 or nopenalty functions. For the `nopenalty' function, the `short.output' should be FALSE.}

\item{y}{data.}

\item{plot_type}{takes two values, 'info.data' or 'info.pair'. 'info.data' is the heatmap of the data with informative variables; 'info.pair' indicates which globally informative variable is pairwise informative for each pair of clusters.}

\item{eps.diff}{The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0. The default value is 1e-5.}

\item{margins}{parameter in 'heatmap.2' function, 2-dimensional numeric vector containing the margins for column and row names, respectively.}

\item{cexRow}{parameter in 'heatmap.2' function, positive numbers for the row axis labeling.}

\item{cexCol}{parameter in 'heatmap.2' function, positive numbers for the column axis labeling.}

\item{lhei}{parameters in 'heatmap.2' function, visual layout of column height.}

\item{lwid}{parameters in 'heatmap.2' function, visual layout of column weight.}

\item{adjCol}{parameters in 'heatmap.2' function, justification of column labels (variables names).}

\item{sepwidth}{parameters in 'heatmap.2' function, 2-dimensional vector giving the width and height of the separator box}
}
\value{
heatmap of the data with informative variables or heatmap of whether the globally informative variables are pairwise informative for each pair of clusters or not.
}
\description{
Heatmaps of the data with estimated informative variables and the indicator for pairwise informativeness of each globally informative variables.
}
\examples{
y <- rbind(matrix(rnorm(120,0,1),ncol=4),
matrix(rnorm(120,4,1), ncol=4), matrix(rnorm(120,0,1),ncol=4))
output <- parse(K = 3, lambda = 1, y=y)
output$mu.hat.best
heatmap_fit(output, y, cexRow=1)

}
\references{
Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2015). gplots: Various R Programming Tools for Plotting Data. R package version 2.17.0.
\url{https://CRAN.R-project.org/package=gplots}
}
\seealso{
\link[gplots]{heatmap.2}
}
\keyword{external}
