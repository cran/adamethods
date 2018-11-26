#' Functional multivariate Frobenius norm
#' 
#' @aliases frobenius_norm_funct_multiv
#'
#' @description 
#' Computes the functional multivariate Frobenius norm.
#' 
#' @usage 
#' frobenius_norm_funct_multiv(m, PM)
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
#' mat <- matrix(1:400, ncol = 20)
#' PM <- matrix(1:100, ncol = 10)
#' frobenius_norm_funct_multiv(mat, PM)
#'                  
#' @export

frobenius_norm_funct_multiv <- function(m, PM){
  di <- dim(m)
  
  s1 <- sum(apply(m[1:(di[1]/2),], 2, int_prod_mat_funct, PM = PM))
  s2 <- sum(apply(m[(di[1]/2 + 1):di[1],], 2, int_prod_mat_funct, PM = PM))
  
  return(s1+s2) 
}