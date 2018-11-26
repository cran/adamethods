#' Archetypoid algorithm with the functional robust Frobenius norm
#' 
#' @aliases archetypoids_funct_robust
#'
#' @description 
#' Archetypoid algorithm with the functional robust Frobenius norm 
#' to be used with functional data. 
#' 
#' @usage 
#' archetypoids_funct_robust(numArchoid, data, huge = 200, ArchObj, PM, prob)
#' 
#' @param numArchoid Number of archetypoids.
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param ArchObj The list object returned by the 
#' \code{\link{stepArchetypesRawData_funct_robust}} function.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
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
#' Irene Epifanio
#' 
#' @seealso 
#' \code{\link[Anthropometry]{archetypoids}}
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919}
#' 
#' @examples 
#' \dontrun{
#' library(fda)
#' ?growth
#' str(growth)
#' hgtm <- t(growth$hgtm)
#' # Create basis:
#' basis_fd <- create.bspline.basis(c(1,ncol(hgtm)), 10)
#' PM <- eval.penalty(basis_fd)
#' # Make fd object:
#' temp_points <- 1:ncol(hgtm)
#' temp_fd <- Data2fd(argvals = temp_points, y = growth$hgtm, basisobj = basis_fd)
#' data_archs <- t(temp_fd$coefs)
#' 
#' lass <- stepArchetypesRawData_funct_robust(data = data_archs, numArch = 3, 
#'                                            numRep = 5, verbose = FALSE, 
#'                                            saveHistory = FALSE, PM, prob = 0.8)
#' 
#' afr <- archetypoids_funct_robust(3, data_archs, huge = 200, ArchObj = lass, PM, 0.8)
#' str(afr)
#' }                                                          
#'                                                      
#' @export

archetypoids_funct_robust <- function(numArchoid, data, huge = 200, ArchObj, PM, prob){
  
  N = dim(data)[1]
  
  ai <- archetypes::bestModel(ArchObj[[1]])
  
  if (is.null(archetypes::parameters(ai))) {
    stop("No archetypes computed")  
  }else{
    ras <- rbind(archetypes::parameters(ai),data)
    dras <- dist(ras, method = "euclidean", diag = FALSE, upper = TRUE, p = 2)
    mdras <- as.matrix(dras)
    diag(mdras) = 1e+11
  }
  
  ini_arch <- sapply(seq(length = numArchoid), nearestToArchetypes, numArchoid, mdras) 
  #print(ini_arch) 
  if (all(ini_arch > numArchoid) == FALSE) {
    k = 1
    neig <- knn(data, archetypes::parameters(ai), 1:N, k = k)
    indices1 <- attr(neig, "nn.index")
    ini_arch <- indices1[,k]
    
    while (any(duplicated(ini_arch))) {
      k = k + 1  
      neig <- knn(data, archetypes::parameters(ai), 1:N, k = k)
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
  rss_ini <- frobenius_norm_funct_robust(resid, PM, prob) / n
  
  res_def <- swap_funct_robust(ini_arch, rss_ini, huge, numArchoid, x_gvv, n, PM, prob)
  
  return(res_def)
}