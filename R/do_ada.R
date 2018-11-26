#' Run the whole classical archetypoid analysis with the Frobenius norm
#' 
#' @aliases do_ada
#'
#' @description 
#' This function executes the entire procedure involved in the archetypoid analysis.
#' Firstly, the initial vector of archetypoids is obtained using the archetypal 
#' algorithm and finally, the optimal vector of archetypoids is returned.
#' 
#' @usage 
#' do_ada(subset, numArchoid, numRep, huge, compare = FALSE,
#'               vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'               outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#'               method = "adjbox")
#' 
#' @param subset Data to obtain archetypes. In ADALARA this is a subset of the 
#' entire data frame.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param compare Boolean argument to compute the robust residual sum of squares 
#' to compare these results with the ones provided by \code{\link{do_ada_robust}}.
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate").
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases: Final vector of archetypoids.
#' \item alphas: Alpha coefficients for the final vector of archetypoids.
#' \item rss: Residual sum of squares corresponding to the final vector of archetypoids.
#' \item rss_rob: If \code{compare=TRUE}, this is the residual sum of squares using
#' the robust Frobenius norm. Otherwise, NULL.
#' \item resid: Vector with the residuals.
#' \item outliers: Outliers.
#' }
#' 
#' @author 
#' Guillermo Vinue, Irene Epifanio
#' 
#' @seealso 
#' \code{\link{stepArchetypesRawData_norm_frob}}, \code{\link{archetypoids_norm_frob}}
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
#' 
#' @examples 
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
#' res_ada <- do_ada(preproc$data, k, numRep, huge, FALSE, method = "adjbox")
#' str(res_ada)     
#' 
#' res_ada1 <- do_ada(preproc$data, k, numRep, huge, FALSE, 
#'                    vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'                    outl_degree = c("outl_strong", "outl_semi_strong", 
#'                                    "outl_moderate"), method = "toler")
#' str(res_ada1) 
#'                  
#' @export

do_ada <- function(subset, numArchoid, numRep, huge, compare = FALSE, 
                   vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
                   outl_degree = c("outl_strong", "outl_semi_strong", 
                                   "outl_moderate"), method = "adjbox") {
  
  lass <- stepArchetypesRawData_norm_frob(data = subset, numArch = numArchoid, 
                                          numRep = numRep, verbose = FALSE, 
                                          saveHistory = FALSE)
  
  ada_subset <- archetypoids_norm_frob(numArchoid, subset, huge = huge, ArchObj = lass) 
  
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
    rss_rob <- frobenius_norm_robust(resid) / n
  }
  
  return(list(cases = k_subset, alphas = alphas_subset, 
              rss = rss_subset, rss_rob = rss_rob, 
              resid = resid_vect, outliers = outl))
}