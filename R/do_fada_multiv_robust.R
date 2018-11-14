#' Run the whole archetypoid analysis with the functional multivariate robust Frobenius norm
#' 
#' @aliases do_fada_multiv_robust
#'
#' @description 
#' This function executes the entire procedure involved in the functional archetypoid 
#' analysis. Firstly, the initial vector of archetypoids is obtained using the 
#' functional archetypal algorithm and finally, the optimal vector of archetypoids is 
#' returned.
#' 
#' @usage 
#' do_fada_multiv_robust(subset, numArchoid, numRep, huge, prob, compare = FALSE, PM,
#'                       method = "adjbox")
#' 
#' @param subset Data to obtain archetypes. In fadalara this is a subset of the 
#' entire data frame.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} times.
#' @param huge Penalization to solve the convex least squares problem, 
#' see \code{\link[Anthropometry]{archetypoids}}.
#' @param prob Probability with values in [0,1].
#' @param compare Boolean argument to compute the non-robust residual sum of squares 
#' to compare these results with the ones provided by \code{\link{do_fada}}.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param method Method to compute the outliers. So far the only option allowed is 
#' 'adjbox' for using adjusted boxplots for skewed distributions. The use of
#' tolerance intervals might also be explored in the future for the multivariate case.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases: Final vector of archetypoids.
#' \item alphas: Alpha coefficients for the final vector of archetypoids.
#' \item rss: Residual sum of squares corresponding to the final vector of archetypoids.
#' \item rss_non_rob: If \code{compare=TRUE}, this is the residual sum of squares using
#' the non-robust Frobenius norm. Otherwise, NULL.
#' \item resid Vector of residuals.
#' \item outliers: Outliers.
#' \item local_rel_imp Matrix with the local (casewise) relative importance 
#' (in percentage) of each variable for the outlier identification. Only for 
#' the multivariate case. It is relative to the outlier observation itself. 
#' The other observations are not considered for computing this importance. 
#' This procedure works because the functional variables are in the same scale, 
#' after standardizing. Otherwise, it couldn't be interpreted like that.
#' \item margi_rel_imp Matrix with the marginal relative importance of each variable 
#' (in percentage) for the outlier identification. Only for the multivariate case. 
#' In this case, the other points are considered, since the value of the outlier 
#' observation is compared with the remaining points.
#' }
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{stepArchetypesRawData_funct_multiv_robust}}, 
#' \code{\link{archetypoids_funct_multiv_robust}}
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
#' set.seed(2018)
#' res_fada <- do_fada_multiv_robust(subset = Xs, numArchoid = 3, numRep = 5, huge = 200, 
#'                                   prob = 0.75, compare = FALSE, PM = PM, method = "adjbox")
#' str(res_fada)
#' res_fada$cases
#' #[1]  8 24 29
#' res_fada$rss
#' #[1] 2.301741
#' }
#'                                   
#' @export

do_fada_multiv_robust <- function(subset, numArchoid, numRep, huge, prob, compare = FALSE, PM, 
                                  method = "adjbox") {
   
  nbasis <- dim(subset)[2] # number of basis.
  nvars <- dim(subset)[3] # number of variables.
  
  lass <- stepArchetypesRawData_funct_multiv_robust(data = subset, numArch = numArchoid, 
                                                    numRep = numRep, verbose = FALSE, 
                                                    saveHistory = FALSE, PM, prob, nbasis, nvars)
  
  fada_subset <- archetypoids_funct_multiv_robust(numArchoid, subset, huge = huge, 
                                                  ArchObj = lass, PM, prob) 
  
  k_subset <- fada_subset$cases # The same with the S3 method anthrCases(fada_subset)
  alphas_subset <- fada_subset$alphas
  
  # Outliers:
  # t(ada_subset$resid) is the residuals matrix with the same dimension
  # as the original data matrix.
  if (method == "adjbox") {
    seq_pts <- sort(c(seq(1, nbasis * nvars, by = nbasis), 
                      rev(nbasis * nvars - nbasis * (1:(nvars-1))), 
                      nbasis * nvars))
    
    odd_pos <- seq(1, length(seq_pts), 2)
    r_list <- list()
    for (i in odd_pos) {
      r_list[[i]] <- apply(fada_subset$resid[seq_pts[i]:seq_pts[i+1],], 2, int_prod_mat_funct, PM = PM)
    }
    r_list1 <- r_list[odd_pos]
    aux <- Reduce(`+`, r_list1)
    resid_vect <- sqrt(aux)
    
    outl_boxb <- boxB(x = resid_vect, k = 1.5, method = method)
    outl <- which(resid_vect > outl_boxb$fences[2]) 
    
    # Local relative importance of the variables in the outlier identification:
    local_rel_imp <- sapply(1:nvars, function(i, j, x, y) round((x[[i]][j]/y[j]) * 100, 2), outl, r_list1, aux)
    if (length(outl) > 0) { # For the case when there are outliers.
      if (length(outl) == 1) { # For the case when there are a single outlier.
        local_rel_imp <- t(local_rel_imp) # We have to create the matrix, otherwise local_rel_imp is a vector.
      }
      dimnames(local_rel_imp) <- list(as.character(outl), paste("V", 1:nvars, sep = ""))
    }
    
    # Marginal relative importance of the variables in the outlier identification:
    margi_rel_imp <- sapply(1:nvars, function(i, j, x) round((x[[i]][j]/sum(x[[i]])) * 100, 2), outl, r_list1)
    if (length(outl) > 0) { # For the case when there are outliers.
      if (length(outl) == 1) { # For the case when there are a single outlier.
        margi_rel_imp <- t(margi_rel_imp) # We have to create the matrix, otherwise margi_rel_imp is a vector.
      }
      dimnames(margi_rel_imp) <- list(as.character(outl), paste("V", 1:nvars, sep = ""))
    }  
  #}else if (method == "toler") {
  #  resid_vect <- apply(fada_subset$resid, 2, int_prod_mat_funct, PM = PM)
  #  # Degree of outlierness:
  #  outl <- do_outl_degree(vect_tol, resid_vect, alpha, 
  #                         paste(outl_degree, "_non_rob", sep = ""))
  }else{
    stop("So far methods allowed are 'adjbox'.") # and 'toler'
  }
  # -----------------------------
  
  rss_subset <- fada_subset$rss
  rss_non_rob <- NULL
  
  if (compare) {
    n <- ncol(t(subset))
    x_gvv <- rbind(t(subset), rep(huge, n))
    
    zs <- x_gvv[,k_subset] 
    zs <- as.matrix(zs)
    
    alphas <- matrix(0, nrow = numArchoid, ncol = n)
    for (j in 1:n) {
      alphas[, j] = coef(nnls(zs, x_gvv[,j]))
    }
    
    resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
    rss_non_rob <- frobenius_norm_funct(resid) / n
  }
  
  #return(c(list(k_subset = k_subset, alphas_subset = alphas_subset, 
  #            rss_subset = rss_subset, rss_non_rob = rss_non_rob), out_tol1))
  return(list(cases = k_subset, alphas = alphas_subset, 
              rss = rss_subset, rss_non_rob = rss_non_rob,
              resid = resid_vect, outliers = outl, 
              local_rel_imp = local_rel_imp, margi_rel_imp = margi_rel_imp))
}