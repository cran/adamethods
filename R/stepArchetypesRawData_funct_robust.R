#' Archetype algorithm to raw data with the functional robust Frobenius norm
#' 
#' @aliases stepArchetypesRawData_funct_robust
#'
#' @description 
#' This is a slight modification of \code{\link[Anthropometry]{stepArchetypesRawData}}
#' to use the functional archetype algorithm with the functional robust Frobenius norm.
#' 
#' @usage 
#' stepArchetypesRawData_funct_robust(data, numArch, numRep = 3, 
#'                                 verbose = TRUE, saveHistory = FALSE, PM, prob)
#' 
#' @param data Data to obtain archetypes.
#' @param numArch Number of archetypes to compute, from 1 to \code{numArch}.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param verbose If TRUE, the progress during execution is shown.
#' @param saveHistory Save execution steps.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
#' 
#' @return 
#' A list with the archetypes.  
#' 
#' @author 
#' Irene Epifanio
#' 
#' @references 
#' Cutler, A. and Breiman, L., Archetypal Analysis. \emph{Technometrics}, 1994,
#' \bold{36(4)}, 338-347, \url{https://doi.org/10.2307/1269949}
#' 
#' Epifanio, I., Functional archetype and archetypoid analysis, 2016. 
#' \emph{Computational Statistics and Data Analysis} \bold{104}, 24-34, 
#' \url{https://doi.org/10.1016/j.csda.2016.06.007}
#' 
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. \emph{Journal of Statistical Software} \bold{30(8)}, 1-23,
#' \url{https://doi.org/10.18637/jss.v030.i08}
#' 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2019. 
#' \emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
#' \url{https://doi.org/10.1016/j.physa.2018.12.036}
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
#' str(lass)   
#' length(lass[[1]])
#' class(lass[[1]])  
#' class(lass[[1]][[5]]) 
#' }                                         
#' 
#' @export

stepArchetypesRawData_funct_robust <- function(data, numArch, numRep = 3, verbose = TRUE, 
                                               saveHistory = FALSE, PM, prob){
  
  mycall <- match.call()
  as <- list()
  for (i in 1:length(numArch)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(numRep)) {
      if (verbose) 
        cat("\n*** numArch=", numArch[i], ", rep=", j, ":\n", sep = "")
        as[[i]][[j]] <- archetypes_funct_robust(data, k = numArch[i], saveHistory = FALSE, 
                                                family = archetypesFamily("original",
                                                                   scalefn = no.scalefn, 
                                                                   rescalefn = no.rescalefn,
                                                                   normfn = frobenius_norm_funct_robust), 
                                                PM = PM, prob = prob)
    }
  }
  return(structure(as, class = "stepArchetypes", call = mycall))
}
