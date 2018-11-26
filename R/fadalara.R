#' Functional parallel archetypoid algorithm for large applications (FADALARA)
#' 
#' @aliases fadalara
#'
#' @description 
#' The FADALARA algorithm is based on the CLARA clustering algorithm. This is the 
#' parallel version of the algorithm. It allows to detect anomalies (outliers). 
#' In the univariate case, there are two different methods to detect them: 
#' the adjusted boxplot (default and most reliable option) and tolerance intervals.
#' In the multivariate case, only adjusted boxplots are used.
#' If needed, tolerance intervals allow to define a degree of outlierness.
#' 
#' @usage 
#' fadalara(data, N, m, numArchoid, numRep, huge, prob, type_alg = "fada", 
#'          compare = FALSE, PM, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'          outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#'          method = "adjbox", multiv, frame)
#' 
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric. The data must have row names
#' so that the algorithm can identify the archetypoids in every sample.
#' @param N Number of samples.
#' @param m Sample size of each sample.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArch}, run the archetype algorithm \code{numRep} 
#' times.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param prob Probability with values in [0,1].
#' @param type_alg String. Options are 'fada' for the non-robust fadalara algorithm, 
#' whereas 'fada_rob' is for the robust fadalara algorithm.
#' @param compare Boolean argument to compute the robust residual sum of squares 
#' if \code{type_alg = "fada"} and the non-robust if \code{type_alg = "fada_rob"}.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate").
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals. 
#' The tolerance intervals are only computed in the univariate case, i.e.,
#' \code{method='toler'} only valid if \code{multiv=FALSE}.
#' @param multiv Multivariate (TRUE) or univariate (FALSE) algorithm.
#' @param frame Boolean value to indicate whether the frame is 
#' computed (Mair et al., 2017) or not. The frame is made up of a subset of
#' extreme points, so the archetypoids are only computed on the frame. 
#' Low frame densities are obtained when only small portions of the data were extreme.
#' However, high frame densities reduce this speed-up.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases Vector of archetypoids.
#' \item rss Optimal residual sum of squares.
#' \item outliers: Outliers.
#' \item alphas: Matrix with the alpha coefficients.
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
#' Guillermo Vinue, Irene Epifanio
#' 
#' @seealso 
#' \code{\link{do_fada}}, \code{\link{do_fada_robust}}
#' 
#' @references 
#' Epifanio, I., Functional archetype and archetypoid analysis, 2016. 
#' \emph{Computational Statistics and Data Analysis} \bold{104}, 24-34, 
#' \url{https://doi.org/10.1016/j.csda.2016.06.007}
#' 
#' Mair, S., Boubekki, A. and Brefeld, U., Frame-based Data Factorizations, 2017.
#' Proceedings of the 34th International Conference on Machine Learning, 
#' Sydney, Australia, 1-9.
#' 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919}
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
#' # Parallel:
#' # Prepare parallelization (including the seed for reproducibility):
#' library(doParallel)
#' no_cores <- detectCores() - 1
#' no_cores 
#' cl <- makeCluster(no_cores)
#' registerDoParallel(cl)
#' clusterSetRNGStream(cl, iseed = 2018)
#' res_fl <- fadalara(data = data_alg, N = N, m = m, numArchoid = k, numRep = numRep, 
#'                    huge = huge, prob = prob, type_alg = "fada_rob", compare = FALSE, 
#'                    PM = PM, method = "adjbox", multiv = TRUE, frame = FALSE) # frame = TRUE
#' stopCluster(cl)
#'
#' res_fl_copy <- res_fl
#' res_fl <- res_fl[which.min(unlist(sapply(res_fl, function(x) x[2])))][[1]]
#' str(res_fl)
#' res_fl$cases
#' #[1]  7  9 10
#' res_fl$rss
#' #[1] 5.772298
#' as.vector(res_fl$outliers)
#' #integer(0)
#' }
#' 
#' @importFrom foreach foreach %dopar%
#'                                                      
#' @export

fadalara <- function(data, N, m, numArchoid, numRep, huge, prob, type_alg = "fada", 
                     compare = FALSE, PM, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
                     outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
                     method = "adjbox", multiv, frame){
  nbasis <- dim(data)[2] # number of basis.
  nvars <- dim(data)[3] # number of variables.
  
  i <- NULL
  n <- nrow(data)
  #rss_aux <- list()
  #k_aux <- list()
  #out_tol <- list()
  rss_aux <- Inf
  rand_obs_iter <- c()
  res <- foreach(i = 1:N, 
          .packages = c("Anthropometry", "archetypes", "nnls", 
                        "adamethods", "tolerance", "utils", "univOutl"))  %dopar% { 
           # Generate subset:
           if (is.null(rand_obs_iter)) {
             #set.seed(seed)
             #set.seed(1)
             rand_obs_si <- sample(1:n, size = m)    
           }else{
             #set.seed(seed)
             #set.seed(1)
             rand_obs_si <- sample(setdiff(1:n, rand_obs_iter), size = m - numArchoid)
             rand_obs_si <- c(rand_obs_si, k_aux)
           }
          # To accumulate the already sampled individuals:
          rand_obs_iter <- c(rand_obs_iter, rand_obs_si)
                    
          # Use FADA on si:
          if (multiv) {
            si <- apply(data, 2:3, function(x) x[rand_obs_si])
          }else{
            si <- data[rand_obs_si,]  
          }
          # Apply the FADA algorithm on si to compute k_si archetypoids:
          if (type_alg == "fada") {
            fada_si <- do_fada(si, numArchoid, numRep, huge, compare, PM, vect_tol, alpha, 
                               outl_degree, method) 
          }else if (type_alg == "fada_rob") { 
            if (multiv) {
              if (frame) {
                g1 <- t(si[,,1])
                G <- dim(si)[3]
                for (i in 2:G) {
                  g12 <- t(si[,,i])  
                  g1 <- rbind(g1, g12) 
                }
                X <- t(g1)
                si_frame <- frame_in_r(X)
                si <- apply(si, 2:3, function(x) x[si_frame])
              }
              fada_si <- do_fada_multiv_robust(si, numArchoid, numRep, huge, prob, compare, PM, method)
            }else{
              if (frame) {
                si_frame <- frame_in_r(si)
                si <- si[si_frame,]
              } 
              fada_si <- do_fada_robust(si, numArchoid, numRep, huge, prob, compare, PM, vect_tol, alpha, 
                                        outl_degree, method) 
            }
          }else{
            stop("Algorithms available are 'fada' or 'fada_rob'.")
          }
                    
           k_si <- fada_si$cases
           alphas_si <- fada_si$alphas
           colnames(alphas_si) <- rownames(si)
           
           # For every observation in data, compute alpha_{k_si} and RSS_{k s_i}:
           if (multiv) {
             rss_si <- do_alphas_rss_multiv(data, si, huge, k_si, rand_obs_si, alphas_si, 
                                            type_alg, PM, prob, nbasis, nvars)
           }else{
             rss_si <- do_alphas_rss(data, si, huge, k_si, rand_obs_si, alphas_si, type_alg, PM, prob)
           }
           
           if (rss_si[[1]] < rss_aux) {
             rss_aux <- rss_si[[1]]
             # Remember to get the right position of the archetypoids regarding the whole data:
             k_aux <- which(rownames(data) %in% rownames(si)[k_si])
             # Matrix with the alpha coefficients:
             alphas_aux <- rss_si[[3]]
             # Residuals:
             resid_aux <- rss_si[[2]]
           }else{
             list(cases = NA, rss = rss_si[[1]], outliers = NA, alphas = NA,
                  local_rel_imp = NA, margi_rel_imp = NA)#,
             #par_samples = rand_obs_si)#, improve = "NO")
           }   
             
           # Compute outliers from the set of archetypoids:
          if (method == "adjbox") {
           if (multiv) {
             seq_pts <- sort(c(seq(1, nbasis*nvars, by = nbasis), 
                               rev(nbasis*nvars - nbasis *(1:(nvars-1))), 
                               nbasis*nvars))
                 
             odd_pos <- seq(1, length(seq_pts), 2)
             r_list <- list()
             for (i in odd_pos) {
               #print(seq_pts[i]:seq_pts[i+1])
               r_list[[i]] <- apply(resid_aux[seq_pts[i]:seq_pts[i+1],], 2, int_prod_mat_funct, PM = PM)
             }
             r_list1 <- r_list[odd_pos]
             aux <- Reduce(`+`, r_list1)
             resid_vect <- sqrt(aux)
           }else{
              resid_vect <- apply(resid_aux, 2, int_prod_mat_sq_funct, PM = PM)
            }
            outl_boxb <- boxB(x = resid_vect, k = 1.5, method = method)
            outl <- which(resid_vect > outl_boxb$fences[2]) 
               
            if (multiv) {
             # Local relative importance of the variables in the outlier identification:
             local_rel_imp <- sapply(1:nvars, function(i, j, x, y) round((x[[i]][j]/y[j]) * 100, 2), 
                                     outl, r_list1, aux)
             if (length(outl) > 0) { # For the case when there are outliers.
              if (length(outl) == 1) { # For the case when there are a single outlier.
               local_rel_imp <- t(local_rel_imp) # We have to create the matrix, otherwise 
                                                 # local_rel_imp is a vector.
              }
              dimnames(local_rel_imp) <- list(as.character(outl), paste("V", 1:nvars, sep = ""))
             }
                 
             # Marginal relative importance of the variables in the outlier identification:
             margi_rel_imp <- sapply(1:nvars, function(i, j, x) round((x[[i]][j]/sum(x[[i]])) * 100, 2), 
                                     outl, r_list1)
             if (length(outl) > 0) { # For the case when there are outliers.
              if (length(outl) == 1) { # For the case when there are a single outlier.
                margi_rel_imp <- t(margi_rel_imp) # We have to create the matrix, otherwise 
                                                  # margi_rel_imp is a vector.
              }
               dimnames(margi_rel_imp) <- list(as.character(outl), paste("V", 1:nvars, sep = ""))
             }  
           }else{
             local_rel_imp <- NULL
             margi_rel_imp <- NULL
           }  
               
         }else if (method == "toler" & multiv == FALSE) {
            resid_vect <- apply(resid_aux, 2, int_prod_mat_funct, PM = PM)
            # Degree of outlierness:
            outl <- do_outl_degree(vect_tol, resid_vect, alpha, 
                                   paste(outl_degree, "_non_rob", sep = ""))
          }
        
           list(cases = k_aux, rss = rss_aux, outliers = outl, alphas = alphas_aux,
                local_rel_imp = local_rel_imp, margi_rel_imp = margi_rel_imp)#,
                  #par_samples = rand_obs_si )#, improve = "YES")

      }
  return(res)
}  