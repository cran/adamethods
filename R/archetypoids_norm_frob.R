#' Archetypoid algorithm with the Frobenius norm
#' 
#' @aliases archetypoids_norm_frob
#'
#' @description 
#' This function is the same as \code{\link[Anthropometry]{archetypoids}} but the 2-norm 
#' is replaced by the Frobenius norm. Thus, the comparison with the robust archetypoids
#' can be directly made.
#' 
#' @usage 
#' archetypoids_norm_frob(numArchoid, data, huge = 200, ArchObj)
#' 
#' @param numArchoid Number of archetypoids.
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param ArchObj The list object returned by the 
#' \code{\link{stepArchetypesRawData_norm_frob}} function. 
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
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. \emph{Journal of Statistical Software} \bold{30(8)}, 1-23,
#' \url{https://doi.org/10.18637/jss.v030.i08}
#' 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2019. 
#' \emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
#' \url{https://doi.org/10.1016/j.physa.2018.12.036}
#' 
#' Vinue, G., Epifanio, I., and Alemany, S., Archetypoids: a new approach to 
#' define representative archetypal data, 2015.
#' \emph{Computational Statistics and Data Analysis} \bold{87}, 102-115,
#' \url{https://doi.org/10.1016/j.csda.2015.01.018}
#' 
#' Vinue, G., Anthropometry: An R Package for Analysis of Anthropometric Data, 2017.
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1-39,
#' \url{https://doi.org/10.18637/jss.v077.i06}
#' 
#' @examples 
#' data(mtcars)
#' data <- mtcars
#' 
#' k <- 3
#' numRep <- 2
#' huge <- 200
#' 
#' lass <- stepArchetypesRawData_norm_frob(data = data, numArch = k, 
#'                                         numRep = numRep, verbose = FALSE)
#' 
#' res <- archetypoids_norm_frob(k, data, huge, ArchObj = lass)
#' str(res)  
#' res$cases
#' res$rss                                                           
#'
#' @importFrom Anthropometry nearestToArchetypes                  
#' @importFrom nnls nnls
#' @importFrom stats coef dist
#' @importFrom FNN knn
#'                                                      
#' @export

archetypoids_norm_frob <- function(numArchoid, data, huge = 200, ArchObj){
  
  N = dim(data)[1]
  
  ai <- archetypes::bestModel(ArchObj[[1]])
  
  # archetypes::parameters(ai) contains the k*m matrix of archetypes.
  if (is.null(archetypes::parameters(ai))) {
   stop("No archetypes computed")  
  }else{
    ras <- rbind(archetypes::parameters(ai),data)
    dras <- dist(ras, method = "euclidean", diag = FALSE, upper = TRUE, p = 2)
    mdras <- as.matrix(dras)
    diag(mdras) = 1e+11
   }
   
  ini_arch <- sapply(seq(length = numArchoid), nearestToArchetypes, numArchoid, mdras) 
    
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
  rss_ini <- frobenius_norm(resid) / n
  #resid <- zs %*% alphas - x_gvv
  #rss_ini <- max(svd(resid)$d) / n
  
  res_def <- swap_norm_frob(ini_arch, rss_ini, huge, numArchoid, x_gvv, n)
  
  return(res_def)
}