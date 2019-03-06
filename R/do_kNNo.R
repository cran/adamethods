#' kNN for outlier detection
#' 
#' @aliases do_knno
#'
#' @description 
#' Ramaswamy et al. proposed the k-nearest neighbors outlier detection method (kNNo). 
#' Each point's anomaly score is the distance to its kth nearest neighbor in the 
#' data set. Then, all points are ranked based on this distance. The higher an example's 
#' score is, the more anomalous it is.
#' 
#' @usage 
#' do_knno(data, k, top_n)
#'
#' @param data Data observations.
#' @param k Number of neighbors of a point that we are interested in.
#' @param top_n Total number of outliers we are interested in.
#' 
#' @return 
#' Vector of outliers.
#' 
#' @author 
#' Guillermo Vinue
#'
#' @references 
#' Ramaswamy, S., Rastogi, R. and Shim, K. 
#' Efficient Algorithms for Mining Outliers from Large Data Sets.
#' SIGMOD'00 Proceedings of the 2000 ACM SIGMOD international conference 
#' on Management of data, 2000, 427-438.  
#'
#' @examples 
#' data(mtcars)
#' data <- as.matrix(mtcars)
#' outl <- do_knno(data, 3, 2)
#' outl
#' data[outl,]
#'                  
#' @export

do_knno <- function(data, k, top_n){
  n <- nrow(data)
  D <- as.matrix(dist(data)) 

  kdist <- c()
  for (i in 1:n) {
    kdist[i] <- (sort(D[i,]))[k + 1]
  }
  anom_score_order <- order(kdist, decreasing = TRUE)
  
  if (top_n != 0) {
    anom_score_order1 <- anom_score_order[1:top_n]
  }else{
    anom_score_order1 <- NULL
  }
  
  return(anom_score_order1)
}