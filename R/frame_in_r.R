#' Compute archetypes frame
#' 
#' @aliases frame_in_r
#'
#' @description 
#' Computing the frame with the approach by Mair et al. (2017).
#' 
#' @usage frame_in_r(X)
#' 
#' @param X Data frame.
#' 
#' @return 
#' Vector with the observations that belong to the frame.
#' 
#' @author 
#' Sebastian Mair, code kindly provided by him.
#' 
#' @references 
#' Mair, S., Boubekki, A. and Brefeld, U., Frame-based Data Factorizations, 2017.
#' Proceedings of the 34th International Conference on Machine Learning, 
#' Sydney, Australia, 1-9.
#' 
#' @examples 
#' \dontrun{
#' X <- mtcars
#' q <- frame_in_r(X)
#' H <- X[q,]
#' q
#' }
#' 
#' @export

frame_in_r <- function(X){
  n <- dim(X)[1]
  Q <- as.matrix( cbind(X,rep(1,n)) )
  ind = c()
  for (i in 1:n) {
    s <- coef(nnls::nnls(t(Q), Q[i,] ))
    ind <- c(ind, which(s != 0))
  }
  ind = sort(unique(ind))
  return(ind)
}
