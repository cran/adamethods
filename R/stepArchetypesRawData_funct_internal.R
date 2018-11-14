#' Perform functional archetypal analysis on a data matrix.
#'
#' @param data A numeric \eqn{n \times m} data matrix.
#' @param k The number of archetypes.
#' @param weights Data weights matrix or vector (used as elements of
#'   the diagonal weights matrix).
#' @param maxIterations The maximum number of iterations.
#' @param minImprovement The minimal value of improvement between two
#'   iterations.
#' @param maxKappa The limit of kappa to report an ill-ness warning.
#' @param verbose Print some details during execution.
#' @param saveHistory Save each execution step in an environment for
#'   further analyses.
#' @param family Blocks defining the underlying problem solving mechanisms;
#'   see \code{\link{archetypesFamily}}.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param ... Additional arguments for family blocks.
#'
#' @return An object of class \code{archetypes}, see
#'   \code{\link{as.archetypes}}.
#'
#'
#' @references 
#' Cutler and Breiman. Archetypal Analysis. Technometrics,
#' 36(4), 1994. 338-348.
#'
#' @examples
#' # Similarly with archetypes_funct.
#' library(archetypes)
#' data(toy)
#' a <- archetypes(toy, 3)
#' str(a)
#'
#' @noRd


archetypes_funct <- function(data, k, weights = NULL, maxIterations = 100, 
                             minImprovement = sqrt(.Machine$double.eps), 
                             maxKappa = 1000, verbose = FALSE, saveHistory = FALSE, 
                             family = archetypesFamily("original"), PM = PM, ...) 
{
  mycall <- match.call()
  famargs <- list(...)
  memento <- NULL
  snapshot <- function(i) {
    a <- list(archetypes = as.archetypes(t(family$rescalefn(x, family$undummyfn(x, zs))), 
                                         k, alphas = t(alphas), betas = t(betas), rss = rss, 
                                         kappas = kappas, 
                                         zas = t(family$rescalefn(x, family$undummyfn(x, zas))), 
                                         residuals = resid, reweights = reweights, weights = weights, 
                                         family = list(class = family$class)))
    memento$save(i, a)
  }
  printIter <- function(i) {
    cat(i, ": rss = ", formatC(rss, 8, format = "f"), ", improvement = ", 
        formatC(imp, 8, format = "f"), "\n", sep = "")
  }
  x1 <- t(data)
  x1 <- family$scalefn(x1, ...)
  x1 <- family$dummyfn(x1, ...)
  x0 <- family$globweightfn(x1, weights, ...)
  x <- x0
  n <- ncol(x)
  m <- nrow(x)
  init <- family$initfn(x, k, ...)
  betas <- init$betas
  alphas <- init$alphas
  zas <- NULL
  zs <- x %*% betas
  #resid <- zs %*% alphas - x
  resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x[1:(nrow(x) - 1),]
  rss <- family$normfn(resid, PM, ...)/n
  reweights <- rep(1, n)
  kappas <- c(alphas = kappa(alphas), betas = kappa(betas), 
              zas = -Inf, zs = kappa(zs))
  isIll <- c(kappas) > maxKappa
  errormsg <- NULL
  if (saveHistory) {
    memento <- new.memento()
    snapshot(0)
  }
  i <- 1
  imp <- +Inf
  tryCatch(while ((i <= maxIterations) & (imp >= minImprovement)) {
    reweights <- family$reweightsfn(resid, reweights, ...)
    x <- family$weightfn(x0, reweights, ...)
    alphas <- family$alphasfn(alphas, zs, x, ...)
    zas <- family$zalphasfn(alphas, x, ...)
    resid1n <- zas[1:(nrow(zas) - 1),] %*% alphas - x[1:(nrow(x) - 1),]
    
    rss1 <- family$normfn(resid1n, PM, ...)/n
    kappas[c("alphas", "zas")] <- c(kappa(alphas), kappa(zas))
    betas <- family$betasfn(betas, x, zas, ...)
    zs <- x %*% betas
    kappas[c("betas", "zs")] <- c(kappa(betas), kappa(zs))
    alphas0 <- family$alphasfn(alphas, zs, x0, ...)
    #resid <- zs %*% alphas0 - x0
    resid <- zs[1:(nrow(zs) - 1),] %*% alphas0 - x0[1:(nrow(x0) - 1),]
    
    rss2 <- family$normfn(resid, PM, ...)/n
    imp <- rss - rss2
    rss <- rss2
    kappas <- c(alphas = kappa(alphas), betas = kappa(betas), 
                zas = kappa(zas), zs = kappa(zs))
    isIll <- isIll & (kappas > maxKappa)
    if (verbose) 
      printIter(i)
    if (saveHistory) 
      snapshot(i)
    i <- i + 1
  }, error = function(e) errormsg <<- e)
  if (!is.null(errormsg)) {
    warning("k=", k, ": ", errormsg)
    return(as.archetypes(NULL, k, NULL, NA, iters = i, call = mycall, 
                         history = history, kappas = kappas))
  }
  if (any(isIll)) 
    warning("k=", k, ": ", paste(names(isIll)[isIll], collapse = ", "), 
            " > maxKappa", sep = "")
  alphas <- family$alphasfn(alphas, zs, x1)
  betas <- family$betasfn(betas, x1, zs)
  zs <- family$undummyfn(x1, zs)
  zs <- family$rescalefn(x1, zs)
  resid <- zs %*% alphas - t(data)
  return(as.archetypes(t(zs), k, t(alphas), rss, iters = (i - 1), call = mycall, history = memento, kappas = kappas, 
                       betas = t(betas), family = family, familyArgs = famargs, 
                       residuals = t(resid), weights = weights, reweights = reweights, 
                       scaling = attr(x1, ".Meta")))
}
