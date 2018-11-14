#' Functional robust Frobenius norm
#' 
#' @aliases frobenius_norm_funct_robust
#'
#' @description 
#' Computes the functional robust Frobenius norm.
#' 
#' @usage 
#' frobenius_norm_funct_robust(m, PM, prob)
#' 
#' @param m Data matrix with the residuals. This matrix has 
#' the same dimensions as the original data matrix.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
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
#' 
#' @examples 
#' library(fda)
#' mat <- matrix(1:9, nrow = 3)
#' fbasis <- create.fourier.basis(rangeval = c(1, 32), nbasis = 3)
#' PM <- eval.penalty(fbasis)
#' frobenius_norm_funct_robust(mat, PM, 0.8)
#'                            
#' @export

frobenius_norm_funct_robust <- function(m, PM, prob) {
  r <- apply(m, 2, int_prod_mat_sq_funct, PM = PM)
  return(sum(bisquare_function(r, prob)))
}