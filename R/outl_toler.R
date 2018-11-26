#' Tolerance outliers
#' 
#' @aliases outl_toler
#'
#' @description 
#' Outliers according to a tolerance interval. This function is used by 
#' the archetypoid algorithms to identify the outliers. See the function
#' \code{nptol.int} in package \code{tolerance}.
#' 
#' @usage 
#' outl_toler(p_tol = 0.95, resid_vect, alpha = 0.05)
#'
#' @param p_tol The proportion of observations to be covered by this 
#' tolerance interval.
#' @param resid_vect Vector of n residuals, where n was the number of rows
#' of the data matrix.
#' @param alpha Significance level.
#' 
#' @return 
#' Vector with the outliers.
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{adalara}}, \code{\link{fadalara}}, \code{\link{do_outl_degree}}
#'  
#' @references 
#' Young, D., tolerance: An R package for estimating tolerance intervals, 2010.
#' \emph{Journal of Statistical Software}, \bold{36(5)}, 1-39, 
#' \url{https://doi.org/10.18637/jss.v036.i05} 
#'  
#' @examples 
#' outl_toler(0.95, 1:100, 0.05)
#'                                                      
#' @export

outl_toler <- function(p_tol = 0.95, resid_vect, alpha = 0.05){
  resid_round <- round(resid_vect, 5)
  tol_interv <- nptol.int(resid_round, alpha = alpha, P = p_tol, side = 1) 
  out_tol <- which(resid_round > tol_interv$`1-sided.upper`)
  out_tol <- as.vector(out_tol)
  return(out_tol)
}
