#' Cleaning outliers
#' 
#' @aliases do_clean
#'
#' @description 
#' Cleaning of the most remarkable outliers. This improves the performance of 
#' the archetypoid algorithm since it is not affected by spurious points.
#' 
#' @usage 
#' do_clean(data, num_pts, range = 1.5, out_perc = 80)
#' 
#' @param data Data frame with (temporal) points in the rows and observations in 
#' the columns.
#' @param num_pts Number of temporal points.
#' @param range Same parameter as in function \code{\link{boxplot}}. 
#' A value of 1.5 is enough to detect amplitude and shift outliers, while a value
#' of 3 is needed to detect isolated outliers.
#' @param out_perc Minimum number of temporal points (in percentage) to consider 
#' the observation as an outlier. Needed when \code{range=1.5}.
#' 
#' @return 
#' Numeric vector with the outliers.
#' 
#' @author 
#' Irene Epifanio
#' 
#' @seealso 
#' \code{\link{boxplot}}
#' 
#' @examples 
#' data(mtcars)
#' data <- mtcars
#' num_pts <- ncol(data)
#' do_clean(t(data), num_pts, 1.5, 80)
#' 
#' @importFrom graphics boxplot
#'                   
#' @export

do_clean <- function(data, num_pts, range = 1.5, out_perc = 80) {
  out_bp <- list()
  for (i in 1:num_pts) {
    b <- boxplot(data[i,], range = range, plot = FALSE)
    if (length(b$out) == 0) { # For cases where there are no outliers.
      # ?boxplot, Value out: the values of any data points which lie beyond the extremes of the whiskers.
      next
    }
    out_bp[[i]] <- which(data[i,] %in% b$out) # By default range=1.5.
  }
  tab_out_bp <- table(unlist(out_bp))
  
  if (range == 1.5) {
    out_min <- (out_perc * num_pts) / 100 
    outl_ada_clean <- as.numeric(names(which(tab_out_bp > out_min)))
  }else if (range == 3){
    outl_ada_clean <- as.numeric(names(which(tab_out_bp >= 1)))
  }
  
  return(outl_ada_clean)
}