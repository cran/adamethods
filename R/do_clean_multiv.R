#' Cleaning multivariate functional outliers
#' 
#' @aliases do_clean_multiv
#'
#' @description 
#' Cleaning of the most remarkable multivariate functional outliers. 
#' This improves the performance of the archetypoid algorithm since it 
#' is not affected by spurious points.
#' 
#' @usage 
#' do_clean_multiv(data, num_pts, range = 1.5, out_perc = 80, nbasis, nvars)
#' 
#' @param data Data frame with (temporal) points in the rows and observations in 
#' the columns.
#' @param num_pts Number of temporal points.
#' @param range Same parameter as in function \code{\link{boxplot}}. 
#' A value of 1.5 is enough to detect amplitude and shift outliers, while a value
#' of 3 is needed to detect isolated outliers.
#' @param out_perc Minimum number of temporal points (in percentage) to consider 
#' the observation as an outlier. Needed when \code{range=1.5}.
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
#' 
#' @return 
#' List with the outliers for each variable.
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{boxplot}}
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
#' x1 <- t(Xs[,,1]) 
#' for (i in 2:nvars) { 
#'  x12 <- t(Xs[,,i]) 
#'  x1 <- rbind(x1, x12) 
#' }
#' data_all <- t(x1) 
#'
#' num_pts <- ncol(data_all) / nvars
#' range <- 3 
#' outl <- do_clean_multiv(t(data_all), num_pts, range, out_perc, nbasis, nvars)
#' outl
#' }
#'                   
#' @export

do_clean_multiv <- function(data, num_pts, range = 1.5, out_perc = 80, nbasis, nvars) {
  seq_pts <- sort(c(seq(1, nbasis*nvars, by = nbasis), 
                    rev(nbasis*nvars - nbasis *(1:(nvars-1))), 
                    nbasis*nvars))
  
  odd_pos <- seq(1, length(seq_pts), 2)
  
  tab_out_bp <- list() ; outl_ada_clean <- list()
  for (i in odd_pos) {
    data_filt <- data[seq_pts[i]:seq_pts[i+1],] 
  
    out_bp <- list()
    for (j in 1:num_pts) {
      b <- boxplot(data_filt[j,], range = range, plot = FALSE)
      if (length(b$out) == 0) { # For cases where there are no outliers.
       # ?boxplot, Value out: the values of any data points which lie beyond the extremes of the whiskers.
       next
      }
      out_bp[[j]] <- which(data_filt[j,] %in% b$out) # By default range=1.5.
    }
    tab_out_bp[[i]] <- table(unlist(out_bp))
  
    if (range == 1.5) {
     out_min <- (out_perc * num_pts) / 100 
     outl_ada_clean[[i]] <- as.numeric(names(which(tab_out_bp[[i]] > out_min)))
    }else if (range == 3){
      outl_ada_clean[[i]] <- as.numeric(names(which(tab_out_bp[[i]] >= 1)))
     }
  }
  outl_ada_clean1 <- outl_ada_clean[odd_pos]
  
  return(outl_ada_clean1)
}  


