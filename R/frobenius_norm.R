#' Frobenius norm
#' 
#' @aliases frobenius_norm
#'
#' @description 
#' Computes the Frobenius norm.
#' 
#' @usage 
#' frobenius_norm(m)
#' 
#' @param m Data matrix with the residuals. This matrix has 
#' the same dimensions as the original data matrix.
#'
#' @details 
#' Residuals are vectors. If there are p variables (columns),
#' for every observation there is a residual that there is 
#' a p-dimensional vector. If there are n observations, the
#' residuals are an n times p matrix. 
#'
#' @return 
#' Real number.
#' 
#' @author 
#' Guillermo Vinue
#
#' @examples 
#' mat <- matrix(1:4, nrow = 2)
#' frobenius_norm(mat)
#'                  
#' @export

frobenius_norm <- function(m) {
  return(sum(apply(m, 2, int_prod_mat)))
}