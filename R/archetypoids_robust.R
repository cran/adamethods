#' Archetypoid algorithm with the robust Frobenius norm
#' 
#' @aliases archetypoids_robust
#'
#' @description 
#' Robust version of the archetypoid algorithm with the Frobenius form.
#' 
#' @usage 
#' archetypoids_robust(numArchoid, data, huge = 200, ArchObj, prob, aaframe)
#' 
#' @param numArchoid Number of archetypoids.
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param ArchObj The list object returned by the 
#' \code{\link{stepArchetypesRawData_robust}} function. 
#' @param prob Probability with values in [0,1].
#' @param aaframe Boolean value to indicate whether the frame-based (TRUE) 
#' (Mair et al., 2017) or the classical (FALSE) (Eugster et al., 2009) archetypes 
#' will be used. The frame-based archetypes are computed with an ancillary python
#' code available at \url{https://www.uv.es/vivigui/software}.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases: Final vector of archetypoids.
#' \item rss: Residual sum of squares corresponding to the final vector of archetypoids.
#' \item archet_ini: Vector of initial archetypoids.
#' \item alphas: Alpha coefficients for the final vector of archetypoids.
#' \item resid: Matrix with the residuals.
#' }
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @importFrom FNN knn
#' 
#' @seealso 
#' \code{\link{archetypoids_norm_frob}}
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
#' data(mtcars)
#' data <- mtcars
#'
#' k <- 3
#' numRep <- 2
#' huge <- 200
#' 
#' lass <- stepArchetypesRawData_robust(data = data, numArch = k, 
#'                                      numRep = numRep, verbose = FALSE, 
#'                                      saveHistory = FALSE, prob = 0.8)
#' 
#' # Please see do_ada_robust.R for an example with aaframe = TRUE.
#' res <- archetypoids_robust(k, data, huge, ArchObj = lass, 0.8, FALSE)
#' str(res)    
#' res$cases
#' res$rss                                                           
#'                  
#' @export

archetypoids_robust <- function(numArchoid, data, huge = 200, ArchObj, prob, aaframe){

  N = dim(data)[1]
  
  if (aaframe) {
    z_frame <- ArchObj
  }else{
    ai <- archetypes::bestModel(ArchObj[[1]])
    z_frame <- archetypes::parameters(ai)
  }
    
  if (is.null(z_frame)) {
   stop("No archetypes computed")  
  }else{
    ras <- rbind(z_frame,data)
    dras <- dist(ras, method = "euclidean", diag = F, upper = T, p = 2)
    mdras <- as.matrix(dras)
    diag(mdras) = 1e+11
  }
    
  ini_arch <- sapply(seq(length = numArchoid), nearestToArchetypes, numArchoid, mdras) 
  
  if (all(ini_arch > numArchoid) == FALSE) {
    k = 1
    neig <- knn(data, z_frame, 1:N, k = k)
    indices1 <- attr(neig, "nn.index")
    ini_arch <- indices1[,k]
    
    while (any(duplicated(ini_arch))) {
      k = k + 1  
      neig <- knn(data, z_frame, 1:N, k = k)
      indicesk <- attr(neig, "nn.index")
      
      dupl <- anyDuplicated(indices1[,1])
      ini_arch <- c(indices1[-dupl,1],indicesk[dupl,k])
    }
  }
  
  n <- ncol(t(data))
  x_gvv <- rbind(t(data), rep(huge, n))
  
  zs <- x_gvv[,ini_arch] 
  zs <- as.matrix(zs)
  
  alphas <- matrix(0, nrow = numArchoid, ncol = n)
  for (j in 1:n) {
    alphas[, j] = coef(nnls(zs, x_gvv[,j]))
  }
  
  resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
  rss_ini <- frobenius_norm_robust(resid, prob) / n
  
  res_def <- swap_robust(ini_arch, rss_ini, huge, numArchoid, x_gvv, n, prob)
  
  return(res_def)
}