#' Archetype algorithm to raw data with the functional multivariate Frobenius norm
#' 
#' @aliases stepArchetypesRawData_funct_multiv
#'
#' @description 
#' This is a slight modification of \code{\link[Anthropometry]{stepArchetypesRawData}}
#' to use the functional archetype algorithm with the multivariate Frobenius norm.
#' 
#' @usage 
#' stepArchetypesRawData_funct_multiv(data, numArch, numRep = 3, 
#'                                    verbose = TRUE, saveHistory = FALSE, PM)
#' 
#' @param data Data to obtain archetypes.
#' @param numArch Number of archetypes to compute, from 1 to \code{numArch}.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param verbose If TRUE, the progress during execution is shown.
#' @param saveHistory Save execution steps.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' 
#' @return 
#' A list with the archetypes. 
#' 
#' @author 
#' Guillermo Vinue
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
#' str(lass)   
#' length(lass[[1]])
#' class(lass[[1]])  
#' class(lass[[1]][[5]])                                             
#' }                                         
#' 
#' @export

stepArchetypesRawData_funct_multiv <- function(data, numArch, numRep = 3, verbose = TRUE, 
                                        saveHistory = FALSE, PM){
  
  mycall <- match.call()
  as <- list()
  for (i in 1:length(numArch)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(numRep)) {
      if (verbose) 
        cat("\n*** numArch=", numArch[i], ", rep=", j, ":\n", sep = "")
        as[[i]][[j]] <- archetypes_funct_multiv(data, k = numArch[i], saveHistory = FALSE, 
                                             family = archetypesFamily("original",
                                                          scalefn = no.scalefn, 
                                                          rescalefn = no.rescalefn,
                                                          normfn = frobenius_norm_funct_multiv,
                                                          dummyfn = make.dummyfn(200)), 
                                             PM = PM)
    }
  }
  return(structure(as, class = "stepArchetypes", call = mycall))
}
