#' Run the whole functional archetypoid analysis with the Frobenius norm
#' 
#' @aliases do_fada
#'
#' @description 
#' This function executes the entire procedure involved in the functional archetypoid 
#' analysis. Firstly, the initial vector of archetypoids is obtained using the 
#' functional archetypal algorithm and finally, the optimal vector of archetypoids is 
#' returned.
#' 
#' @usage 
#' do_fada(subset, numArchoid, numRep, huge, compare = FALSE, PM,
#'               vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'               outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#'               method = "adjbox", prob)
#' 
#' @param subset Data to obtain archetypes. In fadalara this is a subset of the 
#' entire data frame.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param compare Boolean argument to compute the robust residual sum of squares 
#' to compare these results with the ones provided by \code{\link{do_fada_robust}}.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate").
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals.
#' @param prob If \code{compare=TRUE}, probability with values in [0,1].
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases: Final vector of archetypoids.
#' \item alphas: Alpha coefficients for the final vector of archetypoids.
#' \item rss: Residual sum of squares corresponding to the final vector of archetypoids.
#' \item rss_rob: If \code{compare_robust=TRUE}, this is the residual sum of squares using
#' the robust Frobenius norm. Otherwise, NULL.
#' \item resid: Vector of residuals.
#' \item outliers: Outliers.
#' }
#' 
#' @author 
#' Guillermo Vinue, Irene Epifanio
#' 
#' @seealso 
#' \code{\link{stepArchetypesRawData_funct}}, \code{\link{archetypoids_funct}}
#' 
#' @references 
#' Epifanio, I., Functional archetype and archetypoid analysis, 2016. 
#' \emph{Computational Statistics and Data Analysis} \bold{104}, 24-34, 
#' \url{https://doi.org/10.1016/j.csda.2016.06.007}
#' 
#' @examples 
#' \dontrun{
#' library(fda)
#' ?growth
#' str(growth)
#' hgtm <- t(growth$hgtm)
#' 
#' # Create basis:
#' basis_fd <- create.bspline.basis(c(1,ncol(hgtm)), 10)
#' PM <- eval.penalty(basis_fd)
#' # Make fd object:
#' temp_points <- 1:ncol(hgtm)
#' temp_fd <- Data2fd(argvals = temp_points, y = growth$hgtm, basisobj = basis_fd)
#' data_archs <- t(temp_fd$coefs)
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(2018)
#' res_fada <- do_fada(subset = data_archs, numArchoid = 3, numRep = 5, huge = 200, 
#'                     compare = FALSE, PM = PM, method = "adjbox")
#' str(res_fada)      
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(2018)
#' res_fada1 <- do_fada(subset = data_archs, numArchoid = 3, numRep = 5, huge = 200, 
#'                      compare = FALSE, PM = PM, 
#'                      vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'                      outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#'                      method = "toler")
#' str(res_fada1)                               
#' 
#' res_fada2 <- do_fada(subset = data_archs, numArchoid = 3, numRep = 5, huge = 200, 
#'                     compare = TRUE, PM = PM, method = "adjbox", prob = 0.8)
#' str(res_fada2)  
#' 
#' }
#'                                   
#' @export

do_fada <- function(subset, numArchoid, numRep, huge, compare = FALSE, PM, 
                    vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
                    outl_degree = c("outl_strong", "outl_semi_strong", 
                                    "outl_moderate"), method = "adjbox", prob) {
  
  lass <- stepArchetypesRawData_funct(data = subset, numArch = numArchoid, 
                                      numRep = numRep, verbose = FALSE, 
                                      saveHistory = FALSE, PM)
  
  fada_subset <- archetypoids_funct(numArchoid, subset, huge = huge, 
                                    ArchObj = lass, PM) 
  
  k_subset <- fada_subset$cases # The same with the S3 method anthrCases(fada_subset)
  alphas_subset <- fada_subset$alphas
  
  # Outliers:
  # t(ada_subset$resid) is the residuals matrix with the same dimension
  # as the original data matrix.
  if (method == "adjbox") {
    resid_vect <- apply(fada_subset$resid, 2, int_prod_mat_sq_funct, PM = PM)
    outl_boxb <- boxB(x = resid_vect, k = 1.5, method = method)
    outl <- which(resid_vect > outl_boxb$fences[2]) 
  }else if (method == "toler") {
    resid_vect <- apply(fada_subset$resid, 2, int_prod_mat_funct, PM = PM)
    # Degree of outlierness:
    outl <- do_outl_degree(vect_tol, resid_vect, alpha, 
                           paste(outl_degree, "_non_rob", sep = ""))
  }else{
    stop("methods allowed are 'adjbox' and 'toler'.")
  }
  # -----------------------------
  
  rss_subset <- fada_subset$rss
  rss_rob <- NULL
  
  if (compare) {
    n <- ncol(t(subset))
    x_gvv <- rbind(t(subset), rep(huge, n))
    
    zs <- x_gvv[,k_subset] 
    zs <- as.matrix(zs)
    
    alphas <- matrix(0, nrow = numArchoid, ncol = n)
    for (j in 1:n) {
      alphas[, j] = coef(nnls(zs, x_gvv[,j]))
    }
    
    resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
    rss_rob <- frobenius_norm_funct_robust(resid, PM, prob) / n
  }
  
  #return(c(list(k_subset = k_subset, alphas_subset = alphas_subset, 
  #            rss_subset = rss_subset, rss_rob = rss_rob), out_tol1))
  return(list(cases = k_subset, alphas = alphas_subset, 
              rss = rss_subset, rss_rob = rss_rob, 
              resid = resid_vect, outliers = outl))
  
}