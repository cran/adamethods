#' Archetypoid algorithm with the robust Frobenius norm
#' 
#' @aliases archetypoids_robust
#'
#' @description 
#' Robust version of the archetypoid algorithm with the Frobenius form.
#' 
#' @usage 
#' archetypoids_robust(numArchoid, data, huge = 200, ArchObj, prob)
#' 
#' @param numArchoid Number of archetypoids.
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param ArchObj The list object returned by the 
#' \code{\link{stepArchetypesRawData_robust}} function. 
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
#' \code{\link{archetypoids_norm_frob}}
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2019. 
#' \emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
#' \url{https://doi.org/10.1016/j.physa.2018.12.036}
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
#' res <- archetypoids_robust(k, data, huge, ArchObj = lass, 0.8)
#' str(res)    
#' res$cases
#' res$rss                                                           
#'               
#' @importFrom FNN knn                  
#'                        
#' @export

archetypoids_robust <- function(numArchoid, data, huge = 200, ArchObj, prob){

  N = dim(data)[1]
  
  ai <- archetypes::bestModel(ArchObj[[1]])
  z_frame <- archetypes::parameters(ai)
    
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