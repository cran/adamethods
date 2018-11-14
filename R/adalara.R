#' Multivariate parallel archetypoid algorithm for large applications (ADALARA)
#' 
#' @aliases adalara
#'
#' @description 
#' The ADALARA algorithm is based on the CLARA clustering algorithm. This is the 
#' parallel version of the algorithm to try to get faster results. It allows to
#' detect anomalies (outliers). There are two different methods to detect them:
#' the adjusted boxplot (default and most reliable option) and tolerance intervals.
#' If needed, tolerance intervals allow to define a degree of outlierness.
#' 
#' @usage 
#' adalara(data, N, m, numArchoid, numRep, huge, prob, type_alg = "ada", 
#'         compare = FALSE, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#'         outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"), 
#'         method = "adjbox", aaframe, lass)
#' 
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric. The data must have row names
#' so that the algorithm can identify the archetypoids in every sample.
#' @param N Number of samples.
#' @param m Sample size of each sample.
#' @param numArchoid Number of archetypes/archetypoids.
#' @param numRep For each \code{numArchoid}, run the archetype algorithm \code{numRep} 
#' times.
#' @param huge Penalization added to solve the convex least squares problems.
#' @param prob Probability with values in [0,1].
#' @param type_alg String. Options are 'ada' for the non-robust adalara algorithm and 
#' 'ada_rob' for the robust adalara algorithm.
#' @param compare Boolean argument to compute the robust residual sum of squares 
#' if \code{type_alg = "ada"} and the non-robust if \code{type_alg = "ada_rob"}.
#' @param vect_tol Vector with the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate"). 
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals.
#' @param aaframe Boolean value to indicate whether the frame-based (TRUE) 
#' (Mair et al., 2017) or the classical (FALSE) (Eugster et al., 2009) archetypes 
#' will be used. The frame-based archetypes are computed with an ancillary python
#' code available at \url{https://www.uv.es/vivigui/software}.
#' @param lass Frame-based archetypes matrix. Only needed if \code{aaframe = TRUE}.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item cases Optimal vector of archetypoids.
#' \item rss Optimal residual sum of squares.
#' \item outliers: Outliers.
#' }
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{do_ada}}, \code{\link{do_ada_robust}}, \code{\link{adalara_no_paral}}
#' 
#' @references 
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. Journal of Statistical Software 30(8), 1-23.
#' 
#' Mair, S., Boubekki, A. and Brefeld, U., Frame-based Data Factorizations, 2017.
#' Proceedings of the 34th International Conference on Machine Learning, 
#' Sydney, Australia, 1-9.
#' 
#' Vinue, G., (2017). Anthropometry: An R Package for Analysis of Anthropometric Data,
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1--39 
#' 
#' @examples 
#' \dontrun{
#' library(Anthropometry)
#' library(doParallel)
#' 
#' # Prepare parallelization (including the seed for reproducibility):
#' no_cores <- detectCores() - 1
#' cl <- makeCluster(no_cores)
#' registerDoParallel(cl)
#' clusterSetRNGStream(cl, iseed = 1)
#' 
#' # Load data:
#' data(mtcars)
#' data <- mtcars
#' n <- nrow(data)
#' 
#' # Arguments for the archetype/archetypoid algorithm:
#' # Number of archetypoids:
#' k <- 3 
#' numRep <- 2
#' huge <- 200
#' 
#' # Size of the random sample of observations:
#' m <- 10
#' # Number of samples:
#' N <- floor(1 + (n - m)/(m - k))
#' N
#'            
#' prob <- 0.75            
#' 
#' # ADALARA algorithm:
#' preproc <- preprocessing(data, stand = TRUE, percAccomm = 1)
#' data1 <- as.data.frame(preproc$data)
#' 
#' adalara_aux <- adalara(data1, N, m, k, numRep, huge, prob, 
#'                        "ada_rob", FALSE, method = "adjbox", aaframe = FALSE, lass = NULL)
#'
#' #adalara_aux <- adalara(data1, N, m, k, numRep, huge, prob, 
#' #                       "ada_rob", FALSE, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
#' #                       outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
#' #                       method = "toler", aaframe = FALSE, lass = NULL)
#'
#' # Take the minimum RSS, which is in the second position of every sublist:
#' adalara <- adalara_aux[which.min(unlist(sapply(adalara_aux, function(x) x[2])))][[1]]
#' adalara
#'
#' # End parallelization:
#' stopCluster(cl)
#' }
#' 
#' @importFrom foreach foreach %dopar%
#' @importFrom tolerance nptol.int
#' @importFrom utils relist
#'                                                      
#' @export

adalara <- function(data, N, m, numArchoid, numRep, huge, prob, type_alg = "ada", 
                    compare = FALSE, vect_tol = c(0.95, 0.9, 0.85), alpha = 0.05, 
                    outl_degree = c("outl_strong", "outl_semi_strong", "outl_moderate"),
                    method = "adjbox", aaframe, lass){
  
  frame <- NULL
  FurthestSum <- NULL
  ArchetypalAnalysis <- NULL
  
  i <- NULL
  n <- nrow(data)
  #rss_aux <- list()
  #k_aux <- list()
  #out_tol <- list()
  rss_aux <- Inf
  rand_obs_iter <- c()
  res <- foreach(i = 1:N, 
          .packages = c("Anthropometry", "archetypes", "nnls", 
                        "adamethods", "tolerance", "utils"))  %dopar% { # .options.RNG = seed
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
                 
          # Use ADA on si:
          si <- data[rand_obs_si,] 
          # Apply the ADA algorithm on si to compute k_si archetypoids:
          if (type_alg == "ada") {
           ada_si <- do_ada(si, numArchoid, numRep, huge, compare, vect_tol, alpha, 
                            outl_degree, method) 
          }else if (type_alg == "ada_rob") { 
            if (aaframe) {
              X1 <- as.matrix(si)
              X2 <- np_array(X1)      
              q <- frame(X2)  
              
              H <- si[q + 1,] 
              H1 <- as.matrix(H)
              H2 <- np_array(H1)
              Z_init <- H2[FurthestSum(H2, numArchoid)]
              
              PM <- NA
              prob <- NA 
              aa <- ArchetypalAnalysis(H2, Z_init, numArchoid, PM, prob, type_alg, max_iterations = 100, 
                                       stop = FALSE, M = huge, verbose = FALSE, compute_rss = TRUE)
              lass <- aa[[1]] 
              ada_si <- do_ada_robust(H2, numArchoid, numRep, huge, prob, compare, vect_tol, alpha, 
                                      outl_degree, method, aaframe, lass)
            }else{
              ada_si <- do_ada_robust(si, numArchoid, numRep, huge, prob, compare, vect_tol, alpha, 
                                      outl_degree, method, aaframe) 
            }
          }else{
             stop("Algorithms available are 'ada' or 'ada_rob'.")
           }  
                    
           k_si <- ada_si$cases
           alphas_si <- ada_si$alphas
           colnames(alphas_si) <- rownames(si)
           # For every observation in data, compute alpha_{k_si} and RSS_{k s_i}:
           rss_si <- do_alphas_rss(data, si, huge, k_si, rand_obs_si, alphas_si, type_alg, prob = prob)
           
           if (rss_si[[1]] < rss_aux) {
            rss_aux <- rss_si[[1]]
            # Remember to get the right position of the archetypoids regarding the whole data:
            k_aux <- which(rownames(data) %in% rownames(si)[k_si])
                       
            # Compute outliers from the set of archetypoids:
            if (method == "adjbox") {
              resid_vect <- apply(rss_si[[2]], 2, int_prod_mat_sq)
              outl_boxb <- boxB(x = resid_vect, k = 1.5, method = method)
              outl <- which(resid_vect > outl_boxb$fences[2]) 
            }else if (method == "toler") {
              resid_vect <- apply(rss_si[[2]], 2, int_prod_mat)
              # Degree of outlierness:
              outl <- do_outl_degree(vect_tol, resid_vect, alpha, 
                                     paste(outl_degree, "_non_rob", sep = ""))
            }
            
            list(cases = k_aux, rss = rss_aux, outliers = outl)#, improve = "YES")
           }else{
             list(cases = NA, rss = rss_si[[1]], outliers = NA)#, improve = "NO")
            }
           #return(list(cases = k_aux, rss = rss_aux, outliers = out_tol))
         }
  return(res)
}  
