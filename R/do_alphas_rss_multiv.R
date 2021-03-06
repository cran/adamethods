#' Alphas and RSS of every set of multivariate archetypoids
#' 
#' @aliases do_alphas_rss_multiv
#'
#' @description 
#' In the ADALARA algorithm, every time that a set of archetypoids is computed using
#' a sample of the data, the alpha coefficients and the associated residual sum of 
#' squares (RSS) for the entire data set must be computed.
#' 
#' @usage 
#' do_alphas_rss_multiv(data, subset, huge, k_subset, rand_obs, alphas_subset, 
#'                   type_alg = "ada", PM, prob, nbasis, nvars)
#' 
#' @param data Data matrix with all the observations.
#' @param subset Data matrix with a sample of the \code{data} observations.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param k_subset Archetypoids obtained from \code{subset}.
#' @param rand_obs Sample observations that form \code{subset}.
#' @param alphas_subset Alpha coefficients related to \code{k_subset}.
#' @param type_alg String. Options are 'ada' for the non-robust multivariate adalara 
#' algorithm, 'ada_rob' for the robust multivariate adalara algorithm, 'fada' for 
#' the non-robust fda fadalara algorithm and 'fada_rob' for the robust fda
#' fadalara algorithm.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}. Needed when
#' \code{type_alg = 'fada'} or \code{type_alg = 'fada_rob'}.
#' @param prob Probability with values in [0,1]. Needed when
#' \code{type_alg = 'ada_rob'} or \code{type_alg = 'fada_rob'}.
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
#'
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item rss Real number of the residual sum of squares.
#' \item resid_rss Matrix with the residuals.
#' \item alphas Matrix with the alpha values.
#' }
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{archetypoids_norm_frob}}
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
#' # We have to give names to the dimensions to know the 
#' # observations that were identified as archetypoids.
#' dimnames(Xs) <- list(paste("Obs", 1:dim(hgtm)[2], sep = ""), 
#'                      1:nbasis,
#'                      c("boys", "girls"))
#' 
#' n <- dim(Xs)[1] 
#' # Number of archetypoids:
#' k <- 3 
#' numRep <- 20
#' huge <- 200
#'
#' # Size of the random sample of observations:
#' m <- 15
#' # Number of samples:
#' N <- floor(1 + (n - m)/(m - k))
#' N
#' prob <- 0.75
#' data_alg <- Xs
#' 
#' nbasis <- dim(data_alg)[2] # number of basis.
#' nvars <- dim(data_alg)[3] # number of variables.
#' n <- nrow(data_alg)
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(1) 
#' rand_obs_si <- sample(1:n, size = m)  
#' si <- apply(data_alg, 2:3, function(x) x[rand_obs_si])  
#' 
#' fada_si <- do_fada_multiv_robust(si, k, numRep, huge, 0.8, FALSE, PM)
#' 
#' k_si <- fada_si$cases
#' alphas_si <- fada_si$alphas
#' colnames(alphas_si) <- rownames(si)
#' 
#' rss_si <- do_alphas_rss_multiv(data_alg, si, huge, k_si, rand_obs_si, alphas_si, 
#'                                "fada_rob", PM, 0.8, nbasis, nvars)
#' str(rss_si)                                
#' }
#'                                  
#' @export

# Remember that the alphas for the elements belonging to si are already computed.
#subset <- si ; k_subset <- k_si ; rand_obs <- rand_obs_si ; alphas_subset <- alphas_si
do_alphas_rss_multiv <- function(data, subset, huge, k_subset, rand_obs, alphas_subset, 
                                 type_alg = "ada", PM, prob, nbasis, nvars) {
  
  #x11 <- t(data[,,1])
  #x12 <- t(data[,,2])
  #x1 <- rbind(x11,x12)
  x1 <- t(data[,,1])
  B <- dim(data)[3]
  for(i in 2:B){
    x12 <- t(data[,,i])  
    x1 <- rbind(x1, x12) 
  }
  data <- t(x1)
  
  # We have to compute the alphas regarding the archetypoids, so we subset them:
  zs <- t(data)[, which(rownames(data) %in% rownames(subset)[k_subset])] 
  zs <- as.matrix(zs)
  zs <- rbind(zs, huge)
  n <- nrow(data)
  
  # The alphas for the original subset have already been computed, so there is no point
  # in computing them again. We have to compute the alphas only for the non-sampled observations:
  data_remain <- data[setdiff(1:n, rand_obs),]
  ncols <- ncol(t(data_remain))
  data_remain_huge <- rbind(t(data_remain), huge) # The same as rep(huge, ncols)
  
  alphas <- matrix(0, nrow = length(k_subset), ncol = ncols)
  for (j in 1:ncols) {
    alphas[,j] = coef(nnls(zs, data_remain_huge[,j]))
  }
  colnames(alphas) <- rownames(data_remain)
  # Now we join the alphas for the non-sampled and the sampled individuals:
  alphas_all <- rbind(t(alphas), t(alphas_subset))
  # Reorder the rows with the order of the original entire data set, to make it easier
  # the identification of the archetypoids:
  alphas_all1 <- alphas_all[match(rownames(data), rownames(alphas_all)),]
  alphas_all2 <- t(alphas_all1)
  
  # Finally, the rss is computed regarding the entire data set:
  data_huge <- rbind(t(data), huge)
  resid <- zs[1:(nrow(zs) - 1),] %*% alphas_all2 - data_huge[1:(nrow(data_huge) - 1),]
    
  if (type_alg == "ada") {
    rss <- frobenius_norm(resid) / n
  }else if (type_alg == "ada_rob") {
    rss <- frobenius_norm_robust(resid, prob) / n
  }else if (type_alg == "fada") {
    rss <- frobenius_norm_funct(resid, PM) / n
    }else if (type_alg == "fada_rob") {
      rss <- frobenius_norm_funct_multiv_robust(resid, PM, prob, nbasis, nvars) / n
      }else{
    stop("Algorithms available are 'ada', 'ada_rob', 'fada' or 'fada_rob'")
  }  
  
  #resid <- zs %*% alphas_all2 - data_huge
  #rss <- max(svd(resid)$d) / n
  
  #return(rss)
  return(list(rss = rss, resid_rss = resid, alphas = alphas_all2))
}