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
#' Guillermo Vinue, Irene Epifanio
#' 
#' @references 
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. \emph{Journal of Statistical Software} \bold{30(8)}, 1-23,
#' \url{https://doi.org/10.18637/jss.v030.i08}
#' 
#' Vinue, G., Epifanio, I., and Alemany, S.,Archetypoids: a new approach to 
#' define representative archetypal data, 2015.
#' \emph{Computational Statistics and Data Analysis} \bold{87}, 102-115,
#' \url{https://doi.org/10.1016/j.csda.2015.01.018}
#' 
#' Vinue, G., Anthropometry: An R Package for Analysis of Anthropometric Data, 2017.
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1-39,
#' \url{https://doi.org/10.18637/jss.v077.i06}
#
#' @examples 
#' mat <- matrix(1:4, nrow = 2)
#' frobenius_norm(mat)
#'                  
#' @export

frobenius_norm <- function(m) {
  return(sum(apply(m, 2, int_prod_mat)))
}