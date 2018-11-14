#' Degree of outlierness
#' 
#' @aliases do_outl_degree
#'
#' @description 
#' Classification of outliers according to their degree of outlierness.
#' They are classified using the tolerance proportion. For instance,
#' outliers from a 95% tolerance can be considered strong outliers.
#' 
#' @usage 
#' do_outl_degree(vect_tol = c(0.95, 0.9, 0.85), resid_vect, alpha = 0.05, 
#'                outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"))
#'                
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' @param resid_vect Vector of n residuals, where n was the number of rows
#' of the data matrix.
#' @param alpha Significance level. Default 0.05. 
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate"). 
#' 
#' @return 
#' List with the type outliers.
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{outl_toler}}
#'  
#' @examples 
#' do_outl_degree(0.95, 1:100, 0.05, "outl_strong")
#'                                                      
#' @export

do_outl_degree <- function(vect_tol = c(0.95, 0.9, 0.85), resid_vect, alpha = 0.05, 
                           outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate")){
  
  out_tol <- sapply(vect_tol, outl_toler, resid_vect, alpha)
  if (class(out_tol) == "matrix") {
    out_tol <- c(out_tol)
    out_tol1 <- list(out_tol)
  }else{
    # Remove duplicated elements from list:
    un <- unlist(out_tol) 
    out_tol1 <- Map(`[`, out_tol, relist(!duplicated(un), skeleton = out_tol))
  }
  
  names(out_tol1) <- outl_degree
  
  return(out_tol1)
}
