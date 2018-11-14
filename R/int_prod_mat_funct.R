#' Interior product between matrices for FDA
#' 
#' @aliases int_prod_mat_funct
#'
#' @description 
#' Helper function to compute the Frobenius norm 
#' in the functional data analysis (FDA) scenario.
#' 
#' @usage 
#' int_prod_mat_funct(m, PM) 
#' 
#' @param m Data matrix.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' 
#' @return 
#' Data matrix.
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @examples 
#' library(fda)
#' mat <- matrix(1:9, nrow = 3)
#' fbasis <- create.fourier.basis(rangeval = c(1, 32), nbasis = 3)
#' PM <- eval.penalty(fbasis)  
#' int_prod_mat_funct(mat, PM) 
#'                  
#' @export

int_prod_mat_funct <- function(m, PM){
  t(m) %*% PM %*% m
}