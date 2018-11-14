#' Archetypoid algorithm with the functional multivariate Frobenius norm
#' 
#' @aliases archetypoids_funct_multiv
#'
#' @description 
#' Archetypoid algorithm with the functional multivariate Frobenius norm 
#' to be used with functional data. 
#' 
#' @usage 
#' archetypoids_funct_multiv(numArchoid, data, huge = 200, ArchObj, PM)
#' 
#' @param numArchoid Number of archetypoids.
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param ArchObj The list object returned by the 
#' \code{\link{stepArchetypesRawData_funct}} function.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
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
#' @seealso 
#' \code{\link[Anthropometry]{archetypoids}}
#' 
#' @examples 
#' \dontrun{
#' library(fda)
#' ?growth
#' str(growth)
#' hgtm <- growth$hgtm
#' hgtf <- growth$hgtf[,1:39]
#' 
#' # Create array:
#' nvars <- 2
#' data.array <- array(0, dim = c(dim(hgtm), nvars))
#' data.array[,,1] <- as.matrix(hgtm)
#' data.array[,,2] <- as.matrix(hgtf)
#' rownames(data.array) <- 1:nrow(hgtm)
#' colnames(data.array) <- colnames(hgtm)
#' str(data.array)
#' 
#' # Create basis:
#' nbasis <- 10
#' basis_fd <- create.bspline.basis(c(1,nrow(hgtm)), nbasis)
#' PM <- eval.penalty(basis_fd)
#' # Make fd object:
#' temp_points <- 1:nrow(hgtm)
#' temp_fd <- Data2fd(argvals = temp_points, y = data.array, basisobj = basis_fd)
#' 
#' X <- array(0, dim = c(dim(t(temp_fd$coefs[,,1])), nvars))
#' X[,,1] <- t(temp_fd$coef[,,1]) 
#' X[,,2] <- t(temp_fd$coef[,,2])
#' 
#' # Standardize the variables:
#' Xs <- X
#' Xs[,,1] <- scale(X[,,1])
#' Xs[,,2] <- scale(X[,,2])
#' 
#' lass <- stepArchetypesRawData_funct_multiv(data = Xs, numArch = 3, 
#'                                            numRep = 5, verbose = FALSE, 
#'                                            saveHistory = FALSE, PM)
#'                                            
#' afm <- archetypoids_funct_multiv(3, Xs, huge = 200, ArchObj = lass, PM)
#' str(afm)
#' }                                                          
#'                                                      
#' @export

archetypoids_funct_multiv <- function(numArchoid, data, huge = 200, ArchObj, PM){
  
  x11 <- t(data[,,1])
  x12 <- t(data[,,2])
  x1 <- rbind(x11 ,x12)
  data <- t(x1)
  
  N <- dim(data)[1]
  
  ai <- archetypes::bestModel(ArchObj[[1]])
  
  if (is.null(archetypes::parameters(ai))) {
    stop("No archetypes computed")  
  }
  
  dime <- dim(archetypes::parameters(ai))
  ar <- archetypes::parameters(ai)
  R <- matrix(0, nrow = dime[1], ncol = N)
  for (di in 1:dime[1]) { # dime[1] is the number of archetypes.
    Raux <- matrix(ar[di,], byrow = FALSE, ncol = N, nrow = dime[2]) - t(data)
    dii <- dim(Raux)
    R[di,] <- apply(Raux[1:(dii[1]/2),], 2, int_prod_mat_funct, PM = PM) + 
              apply(Raux[(dii[1]/2 + 1):dii[1],], 2, int_prod_mat_funct, PM = PM)
  }
  
  ini_arch <- apply(R, 1, which.min)
  
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
  rss_ini <- frobenius_norm_funct_multiv(resid, PM) / n
  
  res_def <- swap_funct_multiv(ini_arch, rss_ini, huge, numArchoid, x_gvv, n, PM)
  
  return(res_def)
}