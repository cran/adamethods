#' Functional Frobenius norm
#' 
#' @aliases frobenius_norm_funct
#'
#' @description 
#' Computes the functional Frobenius norm.
#' 
#' @usage 
#' frobenius_norm_funct(m, PM)
#' 
#' @param m Data matrix with the residuals. This matrix has 
#' the same dimensions as the original data matrix.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
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
#' Epifanio, I., Functional archetype and archetypoid analysis, 2016. 
#' \emph{Computational Statistics and Data Analysis} \bold{104}, 24-34, 
#' \url{https://doi.org/10.1016/j.csda.2016.06.007}
#' 
#' @examples 
#' library(fda)
#' mat <- matrix(1:9, nrow = 3)
#' fbasis <- create.fourier.basis(rangeval = c(1, 32), nbasis = 3)
#' PM <- eval.penalty(fbasis)
#' frobenius_norm_funct(mat, PM)
#'                  
#' @export

frobenius_norm_funct <- function(m, PM){
  return(sum(apply(m, 2, int_prod_mat_funct, PM = PM))) # 
}