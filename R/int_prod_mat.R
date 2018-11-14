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
#' Guillermo Vinue
#' 
#' @examples 
#' mat <- matrix(1:4, nrow = 2)
#' int_prod_mat(mat)
#' 
#' @export

int_prod_mat <- function(m){
  t(m) %*% m
}