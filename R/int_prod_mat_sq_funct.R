#' Squared interior product between matrices for FDA
#' 
#' @aliases int_prod_mat_sq_funct
#'
#' @description 
#' Helper function to compute the robust Frobenius norm 
#' in the functional data analysis (FDA) scenario.
#' 
#' @usage 
#' int_prod_mat_sq_funct(m, PM) 
#' 
#' @param m Data matrix.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
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
#' library(fda)
#' mat <- matrix(1:9, nrow = 3)
#' fbasis <- create.fourier.basis(rangeval = c(1, 32), nbasis = 3)
#' PM <- eval.penalty(fbasis)  
#' int_prod_mat_sq_funct(mat, PM)
#'                  
#' @export

int_prod_mat_sq_funct <- function(m, PM){
  sqrt(t(m) %*% PM %*% m)
}