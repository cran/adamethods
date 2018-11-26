#' Robust Frobenius norm
#' 
#' @aliases frobenius_norm_robust
#'
#' @description 
#' Computes the robust Frobenius norm.
#' 
#' @usage 
#' frobenius_norm_robust(m, prob)
#' 
#' @param m Data matrix with the residuals. This matrix has 
#' the same dimensions as the original data matrix.
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
#' Irene Epifanio
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919}
#
#' @examples 
#' mat <- matrix(1:4, nrow = 2)
#' frobenius_norm_robust(mat, 0.8)
#'                  
#' @export

frobenius_norm_robust <- function(m, prob){
  r <- apply(m, 2, int_prod_mat_sq)
  return(sum(bisquare_function(r, prob)))
}