#' Multivariate non-parallel archetypoid algorithm for large applications (ADALARA)
#' 
#' @aliases adalara_no_paral
#'
#' @description 
#' The ADALARA algorithm is based on the CLARA clustering algorithm. This is the 
#' non-parallel version of the algorithm. It allows to detect anomalies (outliers). 
#' There are two different methods to detect them: the adjusted boxplot (default 
#' and most reliable option) and tolerance intervals.
#' If needed, tolerance intervals allow to define a degree of outlierness.
#' 
#' @usage 
#' adalara_no_paral(data, seed, N, m, numArchoid, numRep, huge, prob, type_alg = "ada", 
#'                  compare = FALSE, verbose = TRUE, vect_tol = c(0.95, 0.9, 0.85), 
#'                  alpha = 0.05, outl_degree = c("outl_strong", "outl_semi_strong", 
#'                  "outl_moderate"), method = "adjbox", frame)
#' 
#' @param data Data matrix. Each row corresponds to an observation and each column 
#' corresponds to a variable. All variables are numeric. 
#' The data must have row names so that the algorithm can identify the archetypoids 
#' in every sample.
#' @param seed Integer value to set the seed. This ensures reproducibility.
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
#' @param verbose Display progress? Default TRUE.
#' @param vect_tol Vector the tolerance values. Default c(0.95, 0.9, 0.85).
#' Needed if \code{method='toler'}.
#' @param alpha Significance level. Default 0.05. Needed if \code{method='toler'}.
#' @param outl_degree Type of outlier to identify the degree of outlierness.
#' Default c("outl_strong", "outl_semi_strong", "outl_moderate").
#' Needed if \code{method='toler'}.
#' @param method Method to compute the outliers. Options allowed are 'adjbox' for
#' using adjusted boxplots for skewed distributions, and 'toler' for using
#' tolerance intervals.
#' @param frame Boolean value to indicate whether the frame is 
#' computed (Mair et al., 2017) or not. The frame is made up of a subset of
#' extreme points, so the archetypoids are only computed on the frame. 
#' Low frame densities are obtained when only small portions of the data were extreme.
#' However, high frame densities reduce this speed-up.
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
#' Guillermo Vinue, Irene Epifanio
#' 
#' @seealso 
#' \code{\link{do_ada}}, \code{\link{do_ada_robust}}, \code{\link{adalara}}
#'  
#' @references 
#' Eugster, M.J.A. and Leisch, F., From Spider-Man to Hero - Archetypal Analysis in 
#' R, 2009. \emph{Journal of Statistical Software} \bold{30(8)}, 1-23,
#' \url{https://doi.org/10.18637/jss.v030.i08}
#' 
#' Hubert, M. and Vandervieren, E., An adjusted boxplot for skewed distributions, 2008.
#' \emph{Computational Statistics and Data Analysis} \bold{52(12)}, 5186-5201,
#' \url{https://doi.org/10.1016/j.csda.2007.11.008}
#' 
#' Kaufman, L. and Rousseeuw, P.J., Clustering Large Data Sets, 1986.
#' \emph{Pattern Recognition in Practice}, 425-437.
#' 
#' Mair, S., Boubekki, A. and Brefeld, U., Frame-based Data Factorizations, 2017.
#' Proceedings of the 34th International Conference on Machine Learning, 
#' Sydney, Australia, 1-9.
#' 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2019. 
#' \emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
#' \url{https://doi.org/10.1016/j.physa.2018.12.036}
#' 
#' Vinue, G., Anthropometry: An R Package for Analysis of Anthropometric Data, 2017.
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1-39,
#' \url{https://doi.org/10.18637/jss.v077.i06}
#'  
#' @examples 
#' \dontrun{
#' library(Anthropometry)
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
#' res_adalara <- adalara_no_paral(data1, 1, N, m, k, 
#'                                 numRep, huge, prob, "ada_rob", FALSE, TRUE, 
#'                                 method = "adjbox", frame = FALSE)
#'
#' # Examine the results:
#' res_adalara
#' 
#' res_adalara1 <- adalara_no_paral(data1, 1, N, m, k, 
#'                                 numRep, huge, prob, "ada_rob", FALSE, TRUE, 
#'                                 vect_tol = c(0.95, 0.9, 0.85), 
#'                                 alpha = 0.05, outl_degree = c("outl_strong", "outl_semi_strong", 
#'                                                               "outl_moderate"),
#'                                 method = "toler", frame = FALSE)
#' res_adalara1                                 
#'                                  
#' }
#'                                                                                                                                                              
#' @export

adalara_no_paral <- function(data, seed, N, m, numArchoid, numRep, huge, prob, type_alg = "ada", 
                             compare = FALSE, verbose = TRUE, vect_tol = c(0.95, 0.9, 0.85), 
                             alpha = 0.05, outl_degree = c("outl_strong", "outl_semi_strong", 
                                                         "outl_moderate"), 
                             method = "adjbox", frame){

  n <- nrow(data)
  rss_aux <- Inf
  rand_obs_iter <- c()
  for (i in 1:N) {
    if (verbose) {
     print("Iteration:")
     print(i)  
    }
    # Generate subset:
    if (is.null(rand_obs_iter)) {
      set.seed(seed) 
      rand_obs_si <- sample(1:n, size = m)    
    }else{
      set.seed(seed)
      rand_obs_si <- sample(setdiff(1:n, rand_obs_iter), size = m - numArchoid)
      rand_obs_si <- c(rand_obs_si, k_aux)
    }
    # To accumulate the already sampled individuals:
    rand_obs_iter <- c(rand_obs_iter, rand_obs_si)
    #print(rand_obs_si)
    # Use ADA on si:
    si <- data[rand_obs_si,] 
    # Apply the ADA algorithm on si to compute k_si archetypoids:
    if (type_alg == "ada") {
      ada_si <- do_ada(si, numArchoid, numRep, huge, compare, vect_tol, alpha, 
                       outl_degree, method, prob) 
    }else if (type_alg == "ada_rob") { 
      if (frame) {
        si_frame <- frame_in_r(si)
        si <- si[si_frame,]
        rand_obs_si <- rand_obs_si[si_frame]
      }  
      ada_si <- do_ada_robust(si, numArchoid, numRep, huge, prob, compare, vect_tol, alpha, 
                              outl_degree, method)
    }else{
      stop("Algorithms available are 'ada', 'ada_rob'.")
    }
      
    k_si <- ada_si$cases
    alphas_si <- ada_si$alphas
    colnames(alphas_si) <- rownames(si)
    # For every observation in data, compute alpha_{k_si} and RSS_{k s_i}:
    rss_si <- do_alphas_rss(data, si, huge, k_si, rand_obs_si, alphas_si, type_alg, prob = prob)
    
    if (verbose) {
      print("Previous rss value:")
      print(rss_aux)
      print("Current rss value:")
      print(rss_si[[1]])
    }  
    
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
    }  
  } # End loop N iterations.
  return(list(cases = k_aux, rss = rss_aux, outliers = outl))
}
