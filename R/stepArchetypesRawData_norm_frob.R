#' Archetype algorithm to raw data with the Frobenius norm
#' 
#' @aliases stepArchetypesRawData_norm_frob
#'
#' @description 
#' This is a slight modification of \code{\link[Anthropometry]{stepArchetypesRawData}}
#' to use the archetype algorithm with the Frobenius norm.
#' 
#' @usage 
#' stepArchetypesRawData_norm_frob(data, numArch, numRep = 3, 
#'                                 verbose = TRUE, saveHistory = FALSE)
#' 
#' @param data Data to obtain archetypes.
#' @param numArch Number of archetypes to compute, from 1 to \code{numArch}.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param verbose If TRUE, the progress during execution is shown.
#' @param saveHistory Save execution steps.
#' 
#' @return 
#' A list with the archetypes. 
#' 
#' @author 
#' Irene Epifanio
#' 
#' @seealso 
#' \code{\link[Anthropometry]{stepArchetypesRawData}}, 
#' \code{\link[archetypes]{stepArchetypes}}
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
#' data <- as.matrix(mtcars)
#' 
#' numArch <- 5 
#' numRep <- 2
#' 
#' lass <- stepArchetypesRawData_norm_frob(data = data, numArch = 1:numArch, 
#'                                         numRep = numRep, verbose = FALSE)
#'                                         
#' str(lass)   
#' length(lass[[1]])
#' class(lass[[1]])  
#'                                                    
#' @importFrom archetypes archetypes archetypesFamily
#' @importFrom utils history
#' 
#' @export

stepArchetypesRawData_norm_frob <- function(data, numArch, numRep = 3, 
                                            verbose = TRUE, saveHistory = FALSE){
  
  mycall <- match.call()
  as <- list()
  for (i in 1:length(numArch)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(numRep)) {
      if (verbose) 
       cat("\n*** numArch=", numArch[i], ", rep=", j, ":\n", sep = "")
       as[[i]][[j]] <- archetypes_norm_frob(data, k = numArch[i], saveHistory = FALSE,  
                                            family = archetypesFamily("original",
                                                              scalefn = no.scalefn,
                                                              rescalefn = no.rescalefn,
                                                              normfn = frobenius_norm))
    }
  }
  return(structure(as, class = 'stepArchetypes', call = mycall))
}