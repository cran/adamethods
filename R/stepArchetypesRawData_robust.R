#' Archetype algorithm to raw data with the robust Frobenius norm
#' 
#' @aliases stepArchetypesRawData_robust
#'
#' @description 
#' This is a slight modification of \code{\link[Anthropometry]{stepArchetypesRawData}}
#' to use the archetype algorithm with the robust Frobenius norm.
#' 
#' @usage 
#' stepArchetypesRawData_robust(data, numArch, numRep = 3, 
#'                              verbose = TRUE, saveHistory = FALSE, prob)
#' 
#' @param data Data to obtain archetypes.
#' @param numArch Number of archetypes to compute, from 1 to \code{numArch}.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param verbose If TRUE, the progress during execution is shown.
#' @param saveHistory Save execution steps.
#' @param prob Probability with values in [0,1].
#'
#' @return 
#' A list with the archetypes.
#'
#' @author 
#' Irene Epifanio
#' 
#' @seealso 
#' \code{\link{stepArchetypesRawData_norm_frob}}
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919}
#' 
#' @examples 
#' data(mtcars)
#' data <- as.matrix(mtcars)
#' 
#' numArch <- 5 
#' numRep <- 2
#' 
#' lass <- stepArchetypesRawData_robust(data = data, numArch = 1:numArch, 
#'                                      numRep = numRep, verbose = FALSE,
#'                                      saveHistory = FALSE, prob = 0.8)
#' str(lass)   
#' length(lass[[1]])
#' class(lass[[1]])                                       
#'                  
#' @export

stepArchetypesRawData_robust <- function(data, numArch, numRep = 3, verbose = TRUE, 
                                         saveHistory = FALSE, prob){
  
  mycall <- match.call()
  as <- list()
  for (i in 1:length(numArch)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(numRep)) {
      if (verbose) 
       cat("\n*** numArch=", numArch[i], ", rep=", j, ":\n", sep = "")
       as[[i]][[j]] <- archetypes_robust(data, k = numArch[i], saveHistory = FALSE, 
                                         family = archetypesFamily("original",
                                                         scalefn = no.scalefn, 
                                                         rescalefn = no.rescalefn,
                                                         normfn = frobenius_norm_robust),
                                         prob = prob)
    }
  }
  return(structure(as, class = "stepArchetypes", call = mycall))
}