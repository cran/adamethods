#' Bisquare function
#' 
#' @aliases  bisquare_function
#'
#' @description 
#' This function belongs to the bisquare family of loss functions.
#' The bisquare family can better cope with extreme outliers.
#' 
#' @usage 
#' bisquare_function(resid, prob, ...)
#' 
#' @param resid Vector of residuals, computed from the 
#' \eqn{m \times n} residuals data matrix.
#' @param prob Probability with values in [0,1].
#' @param ... Additional possible arguments.
#'  
#' @return 
#' Vector of real numbers.
#'
#' @author 
#' Irene Epifanio
#'
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919} 
#'
#' @examples 
#' resid <- c(2.47, 11.85)  
#' bisquare_function(resid, 0.8)
#' 
#' @importFrom stats median quantile
#'                  
#' @export

bisquare_function <- function(resid, prob, ...) {
  resid0 <- resid < sqrt(.Machine$double.eps)
  c <- quantile(resid[!resid0], probs = prob)
  v <- resid / c
  co <- c^2/6
  
  ifelse(resid <= c, co*(1 - (1 - v^2)^3), co)
}
