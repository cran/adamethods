#' Interior product between matrices
#' 
#' @aliases int_prod_mat
#'
#' @description 
#' Helper function to compute the Frobenius norm.
#' 
#' @usage 
#' int_prod_mat(m) 
#' 
#' @param m Data matrix.
#' 
#' @return 
#' Data matrix.
#' 
#' @author 
#' Irene Epifanio
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2019. 
#' \emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
#' \url{https://doi.org/10.1016/j.physa.2018.12.036}
#' 
#' @examples 
#' mat <- matrix(1:4, nrow = 2)
#' int_prod_mat(mat)
#' 
#' @export

int_prod_mat <- function(m){
  t(m) %*% m
}