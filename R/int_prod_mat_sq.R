#' Squared interior product between matrices
#' 
#' @aliases int_prod_mat_sq
#'
#' @description 
#' Helper function to compute the robust Frobenius norm.
#' 
#' @usage 
#' int_prod_mat_sq(m) 
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
#' int_prod_mat_sq(mat)
#' 
#' @export

int_prod_mat_sq <- function(m){
  sqrt(t(m) %*% m)
}