#' Run the whole robust archetypoid analysis with the robust Frobenius norm
#' 
#' @aliases do_ada_robust
#'
#' @description 
#' This function executes the entire procedure involved in the robust archetypoid 
#' analysis. Firstly, the initial vector of archetypoids is obtained using the 
#' robust archetypal algorithm and finally, the optimal vector of robust archetypoids 
#' is returned.
#' 
#' @usage 
#' do_ada_robust(subset, numArchoid, numRep, huge, prob, compare = FALSE,
#'               vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'               outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#'               method = "adjbox", aaframe, lass)
#' 
#' @param subset Data to obtain archetypes. In ADALARA this is a subset of the 
#' entire data frame.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} 
#' times.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param prob Probability with values in [0,1].
#' @param compare Boolean argument to compute the non-robust residual sum of squares 
#' to compare these results with the ones provided by \code{\link{do_ada}}.
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate").
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals.
#' @param aaframe Boolean value to indicate whether the frame-based (TRUE) 
#' (Mair et al., 2017) or the classical (FALSE) (Eugster et al., 2009) archetypes 
#' will be used. The frame-based archetypes are computed with an ancillary python
#' code available at \url{https://www.uv.es/vivigui/software}.
#' @param lass Frame-based archetypes matrix. Needed if \code{aaframe = TRUE}.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases: Final vector of archetypoids.
#' \item alphas: Alpha coefficients for the final vector of archetypoids.
#' \item rss: Residual sum of squares corresponding to the final vector of archetypoids.
#' \item rss_non_rob: If \code{compare=TRUE}, this is the residual sum of squares using
#' the non-robust Frobenius norm. Otherwise, NULL.
#' \item resid Vector of residuals.
#' \item outliers: Outliers.
#' }
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{stepArchetypesRawData_robust}}, \code{\link{archetypoids_robust}}
#' 
#' @references 
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. Journal of Statistical Software 30(8), 1-23.
#' 
#' Mair, S., Boubekki, A. and Brefeld, U., Frame-based Data Factorizations, 2017.
#' Proceedings of the 34th International Conference on Machine Learning, 
#' Sydney, Australia, 1-9.
#' 
#' Vinue, G., (2017). Anthropometry: An R Package for Analysis of Anthropometric Data,
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1--39  
#' 
#' @examples 
#' \dontrun{
#' library(Anthropometry)
#' data(mtcars)
#' #data <- as.matrix(mtcars)
#' data <- mtcars
#' 
#' k <- 3
#' numRep <- 2
#' huge <- 200
#' 
#' preproc <- preprocessing(data, stand = TRUE, percAccomm = 1)
#' set.seed(2018)
#' res_ada_rob <- do_ada_robust(preproc$data, k, numRep, huge, 0.8,
#'                              FALSE, method = "adjbox", aaframe = FALSE)
#' str(res_ada_rob)    
#' 
#' res_ada_rob1 <- do_ada_robust(preproc$data, k, numRep, huge, 0.8,
#'                              FALSE, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'                              outl_degree = c("outl_strong", "outl_semi_strong", 
#'                                              "outl_moderate"),
#'                              method = "toler", aaframe = FALSE)
#' str(res_ada_rob1)  
#'  
#' # ---------------                          
#'                                                                          
#' library(reticulate)
#' source_python("frame_aa.py") # \url{https://www.uv.es/vivigui/software}.
#' X <- read.csv("USAFSurvey.csv", header = FALSE) # \url{https://www.uv.es/vivigui/software}.
#' X1 <- as.matrix(X)
#' X2 <- np_array(X1)      
#' q <- frame(X2)                 
#' cat(paste("The frame density of USAFSurvey is ", length(q) / nrow(X) * 100, "%", sep = ""))
#' 
#' H <- X[q + 1,] # q+1 to get the right indexes, because python starts counting from 0.
#' H1 <- as.matrix(H)
#' H2 <- np_array(H1)
#' Z_init <- H2[FurthestSum(H2, k)]
#' 
#' PM <- NA
#' prob <- NA 
#' alg <- "ada"
#' aa <- ArchetypalAnalysis(H2, Z_init, k, PM, prob, alg, max_iterations = 100, stop = FALSE, 
#'                          M = huge, verbose = FALSE, compute_rss = TRUE)
#' lass <- aa[[1]]  
#' 
#' res_ada1 <- do_ada_robust(subset = H, numArchoid = k, numRep = numRep, 
#'                           huge = huge, prob =0.8, compare = FALSE, 
#'                           method = "adjbox", aaframe = TRUE, lass = lass)  
#' # It takes a few seconds.                           
#' str(res_ada1)                                                                        
#' }
#'                  
#' @importFrom univOutl boxB
#' 
#' @export

do_ada_robust <- function(subset, numArchoid, numRep, huge, prob, compare = FALSE, 
                          vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
                          outl_degree = c("outl_strong", "outl_semi_strong", 
                                          "outl_moderate"), method = "adjbox", aaframe, lass) {
  
  if (!aaframe) {
    lass <- stepArchetypesRawData_robust(data = subset, numArch = numArchoid, 
                                         numRep = numRep, verbose = FALSE, 
                                         saveHistory = FALSE, prob)
  }  
  
  ada_subset <- archetypoids_robust(numArchoid, subset, huge = huge, ArchObj = lass, prob, aaframe) 
  
  k_subset <- ada_subset$cases # The same with the S3 method anthrCases(ada_subset)
  alphas_subset <- ada_subset$alphas
  
  # Outliers:
  # t(ada_subset$resid) is the residuals matrix with the same dimension
  # as the original data matrix.
  if (method == "adjbox") {
    resid_vect <- apply(ada_subset$resid, 2, int_prod_mat_sq)
    outl_boxb <- boxB(x = resid_vect, k = 1.5, method = method)
    outl <- which(resid_vect > outl_boxb$fences[2]) 
  }else if (method == "toler") {
    resid_vect <- apply(ada_subset$resid, 2, int_prod_mat)
    # Degree of outlierness:
    outl <- do_outl_degree(vect_tol, resid_vect, alpha, 
                           paste(outl_degree, "_non_rob", sep = ""))
  }else{
     stop("methods allowed are 'adjbox' and 'toler'.")
   }
  # -----------------------------
  
  rss_subset <- ada_subset$rss
  rss_non_rob <- NULL
  
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
    rss_non_rob <- frobenius_norm(resid) / n
  }
  
  return(list(cases = k_subset, alphas = alphas_subset, 
              rss = rss_subset, rss_non_rob = rss_non_rob,
              resid = resid_vect, outliers = outl))
}
