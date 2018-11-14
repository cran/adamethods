#' Functional multivariate robust Frobenius norm
#' 
#' @aliases frobenius_norm_funct_multiv_robust
#'
#' @description 
#' Computes the functional multivariate robust Frobenius norm.
#' 
#' @usage 
#' frobenius_norm_funct_multiv_robust(m, PM, prob, nbasis, nvars)
#' 
#' @param m Data matrix with the residuals. This matrix has 
#' the same dimensions as the original data matrix.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
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
#' mat <- matrix(1:400, ncol = 20)
#' PM <- matrix(1:100, ncol = 10)
#' frobenius_norm_funct_multiv_robust(mat, PM, 0.8, 10, 2)
#'                  
#' @export

frobenius_norm_funct_multiv_robust <- function(m, PM, prob, nbasis, nvars){
  seq_pts <- sort(c(seq(1, nbasis*nvars, by = nbasis), 
                    rev(nbasis*nvars - nbasis *(1:(nvars-1))), 
                    nbasis*nvars))
  
  #di <- dim(m)
  #r1 <- apply(m[1:(di[1]/2),], 2, int_prod_mat_funct, PM = PM)
  #r2 <- apply(m[(di[1]/2 + 1):di[1],], 2, int_prod_mat_funct, PM = PM)
  #r <- sqrt(r1+r2)
  
  odd_pos <- seq(1, length(seq_pts), 2)
  r_list <- list()
  for (i in odd_pos) {
    #print(seq_pts[i]:seq_pts[i+1])
    r_list[[i]] <- apply(m[seq_pts[i]:seq_pts[i+1],], 2, int_prod_mat_funct, PM = PM)
  }
  r_list1 <- r_list[odd_pos]
  r <- sqrt(Reduce(`+`, r_list1))
  
  return(sum(bisquare_function(r, prob)))
}