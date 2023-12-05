#' @name response2drug
#' @aliases response2drug
#'
#' @title Gene-expression Data for Asthma Disease
#' @description The data contains gene expression data of 108 objects.
#' @usage data(response2drug)
#' @docType data
#' @format
#' \describe{
#'  \item{Rows}{Objects}
#'  \item{Columns}{Genes}
#' }
#'
#' @details
#' This is the microarray gene expression data from NCBI's Gene Expression Omnibus database with Gene Expression Omnibus Series accession number GSE43696. The data consist of 108 samples with 20 healthy, 50 moderate asthma and 38 severe asthma patients; and 405 genes. Each row represents one observation and each column represents one gene.
#'
#' @source
#' \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43696}
#'
#' @references
#' Voraphani, N., Gladwin, M.T., Contreras, A.U., Kaminski, N., Tedrow, J.R., Milosevic, J., Bleecker, E.R., Meyers, D.A., Ray, A., Ray, P. and Erzurum, S.C. (2014) An airway epithelial iNOS-DUOX2-thyroid peroxidase metabolome drives Th1/Th2 nitrative stress in human severe asthma. \emph{Mucosal immunology} \bold{7(5)}, 1175-1185.
#'
#' @keywords
#' dataset
#' @examples
#' data(response2drug)
#' output1 = parse (K=2, lambda = 1, y = response2drug, N = 100,
#' kms.iter = 50, kms.nstart = 50, eps.diff = 1e-5, eps.em = 1e-5,
#' model.crit = 'gic', backward = TRUE, cores = 2)
#' output1$mu.hat.best[, 1:5]

"response2drug"
