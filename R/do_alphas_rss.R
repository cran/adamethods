#' Alphas and RSS of every set of archetypoids
#' 
#' @aliases do_alphas_rss
#'
#' @description 
#' In the ADALARA algorithm, every time that a set of archetypoids is computed using
#' a sample of the data, the alpha coefficients and the associated residual sum of 
#' squares (RSS) for the entire data set must be computed.
#' 
#' @usage 
#' do_alphas_rss(data, subset, huge, k_subset, rand_obs, alphas_subset, 
#'               type_alg = "ada", PM, prob)
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
#' data(mtcars)
#' data <- mtcars
#' n <- nrow(data)
#' m <- 10
#' 
#' k <- 3 
#' numRep <- 2
#' huge <- 200
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(1)
#' rand_obs_si <- sample(1:n, size = m) 
#' 
#' si <- data[rand_obs_si,]
#' ada_si <- do_ada(si, k, numRep, huge, FALSE) 
#'
#' k_si <- ada_si$cases
#' alphas_si <- ada_si$alphas
#' colnames(alphas_si) <- rownames(si)     
#' 
#' rss_si <- do_alphas_rss(data, si, huge, k_si, rand_obs_si, alphas_si, "ada")
#' str(rss_si)
#'                                  
#' @export

# Remember that the alphas for the elements belonging to si are already computed.
#subset <- si ; k_subset <- k_si ; rand_obs <- rand_obs_si ; alphas_subset <- alphas_si
do_alphas_rss <- function(data, subset, huge, k_subset, rand_obs, alphas_subset, type_alg = "ada", PM, prob) {
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
      rss <- frobenius_norm_funct_robust(resid, PM, prob) / n
      }else{
    stop("Algorithms available are 'ada', 'ada_rob', 'fada' or 'fada_rob'")
  }  
  
  #resid <- zs %*% alphas_all2 - data_huge
  #rss <- max(svd(resid)$d) / n
  
  #return(rss)
  return(list(rss = rss, resid_rss = resid, alphas = alphas_all2))
}