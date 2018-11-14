#' archetypoids_norm_frob_internal
#' 
#' @rdname archetypoids_norm_frob_internal
#'
#' @description 
#' This file contains six extra functions. The \code{ahull} function is a helper 
#' function to calculate the approximated convex hull. The \code{no.scalefn} and 
#' \code{no.rescalefn} functions allow us to use the archetype algorithm without 
#' standardizing the data (see \code{\link{stepArchetypesRawData_norm_frob}}). 
#' The other three functions are called by the archetypoid algorithm 
#' (\code{\link{archetypoids_norm_frob}}). All these functions are therefore not 
#' solely used.
#' 
#' @references 
#' Vinue, G., (2017). Anthropometry: An R Package for Analysis of Anthropometric Data,
#' \emph{Journal of Statistical Software} \bold{77(6)}, 1--39 
#' 
#' Vinue, G., Epifanio, I., and Alemany, S., (2015). Archetypoids: a new approach to define 
#' representative archetypal data, \emph{Computational Statistics and Data Analysis} \bold{87}, 
#' 102--115.
#' 
#' Eugster, M. J., and Leisch, F., (2009). From Spider-Man to Hero - Archetypal Analysis in R, 
#' \emph{Journal of Statistical Software} \bold{30}, 1--23, \url{http://www.jstatsoft.org/}.
#' 
#' @author 
#' Guillermo Vinue
#' 
#' @seealso 
#' \code{\link{archetypoids_norm_frob}}
#'                  

#' Helper function to calculate the approximated convex hull.
#' @param zs An \code{archetypes} object.
#' @return Matrix with the points.
#' @noRd
ahull <- function(zs) {
  a <- rbind(archetypes::parameters(zs), archetypes::parameters(zs)[1,])
  xc <- a[,1]; xm <- mean(xc)
  yc <- a[,2]; ym <- mean(yc)
  
  real <- xc - xm
  imag <- yc - ym
  angle <- atan2(imag, real)
  
  index <- order(angle)
  
  return(a[c(index, index[1]),])
}

#' Function to avoid scaling
#' @param x Data matrix.
#' @return Same data matrix.
#' @noRd
no.scalefn <- function(x, ...){
  return(x)
}

#' Function to avoid rescaling
#' @param x Data matrix.
#' @param zs Archetypes matrix.
#' @return Same data matrix.
#' @noRd
no.rescalefn <- function(x, zs, ...){
  if ( is.null(zs) )
    return(matrix(NA, nrow = 0, ncol = 0))
  return(zs)
}

#' General swap function of the archetypoid analysis 
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss_arch_ini Residual sum of squares associated to \code{vect_arch_ini}.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @return Optimal vector of archetypoids.
#' @noRd
swap_norm_frob <- function(vect_arch_ini, rss_arch_ini, huge = 200, numArchoid, x_gvv, n){
  
  vect_arch_end <- vect_arch_ini
  rss <- rss_arch_ini
  
  
  for (l in 1:numArchoid) {
    
    rss1 <- c()
    
    setpossibles <- setdiff(1:n,vect_arch_ini)
    
    for (i in setpossibles) {
      zs <- x_gvv[,c(i,vect_arch_ini[-l])] 
      zs <- as.matrix(zs)
      alphas <- matrix(0, nrow = numArchoid, ncol = n)
      for (j in 1:n) {
        alphas[, j] = coef(nnls(zs, x_gvv[,j]))
      }
      
      resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
      rss1[i] <- frobenius_norm(resid) / n
      #resid <- zs %*% alphas - x_gvv
      #rss1[i] <- max(svd(resid)$d) / n
      
      if (rss1[i] < rss) {
        rss <- rss1[i]
        vect_arch_end = c(i,vect_arch_ini[-l])
      }
    }
  }
  
  if (numArchoid == 1) {
    result <- swap2_k1_norm_frob(vect_arch_end, vect_arch_ini, rss, huge, numArchoid, x_gvv, n)
  }else{
    result <- swap2_norm_frob(vect_arch_end, vect_arch_ini, rss, huge, numArchoid, x_gvv, n)
  } 
  
  return(result)
}

#' Further swap function of archetypoid analysis 
#' @param vect_arch_end Current vector of archetypoids.
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss Current residual sum of squares of the vector of archetypoids.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @return Optimal vector of archetypoids.
#' @noRd
swap2_norm_frob <- function(vect_arch_end, vect_arch_ini, rss, huge = 200, numArchoid, x_gvv, n){
  
  vect_arch_ini_aux <- vect_arch_ini
  vect_arch_end_aux <- vect_arch_end
  
  rss_aux <- rss
  vect_arch_end2 <- c()
  
  while (any(sort(vect_arch_end_aux) != sort(vect_arch_ini_aux))) { #this loop is executed while both vectors are 
    #different in at least one element. Since we have not found an R function that says whether two vectors are 
    #equal, we have done the following: first, we have ordered them (because both vectors may be the same, although 
    #their elements are in a different order, for example, c(1,2,3) and c(3,1,2)). Then, we have checked if any 
    #element does not match with the element that is placed in the same position in the other vector.
    
    se <- setdiff(vect_arch_end_aux,vect_arch_ini_aux) #this function looks for the distinct element between the 
    #initial vector (in the first iteration, the initial vector is either the nearest or which vector, while in the 
    #second iteration, the initial vector is the final vector returned by the swap step and so on and so forth) 
    #and the final vector (the final vector in the first iteration is that one returned by the swap step but in the 
    #following iterations, it is the vector returned by the swap2 function).
    
    se1 <- setdiff(vect_arch_end_aux,se) #the elements different from the distinct one of the former setdiff 
    #function.
    
    for (l in 1:length(se1)) {
      
      rss1 <- c()
      
      comp <- c(se,se1[-l]) #vector made up of the distinct element with the no distincts without one.
      setpossibles <- setdiff(1:n,vect_arch_end_aux)
      
      for (i in setpossibles) {
        zs <- x_gvv[,c(i,comp)] 
        zs <- as.matrix(zs)
        alphas <- matrix(0, nrow = numArchoid, ncol = n)
        for (j in 1:n) {
          alphas[, j] = coef(nnls(zs, x_gvv[,j]))
        }
        
        resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
        rss1[i] <- frobenius_norm(resid) / n
        #resid <- zs %*% alphas - x_gvv
        #rss1[i] <- max(svd(resid)$d) / n
        
        if (rss1[i] < rss_aux) {
          rss_aux <- rss1[i]
          vect_arch_end2 <- c(i,comp)
        }
      }
    }
    
    if (is.null(vect_arch_end2)) { #if vect_arch_end2 is NULL, this means that any vector improves the final vector 
      #of the swap function (called vect_arch_end_aux). Therefore, the vect_arch_end_aux that it is going to be 
      #returned is just the vect_arch_end_aux that it is the final vector of the first swap. If we don't add this 
      #condition, the function displays an error because it might happen that the vect_arch_end returned by the 
      #swap function already is the best vector (this happens with the nba2d database) and therefore the swap2 
      #function is not going to be able to improve it. 
      vect_arch_end_aux <- vect_arch_end_aux
      vect_arch_ini_aux <- vect_arch_end_aux #In addition, we also have to fix the initial vector as the final 
      #vector in order to the while loop stops.
    }else{
      vect_arch_ini_aux <- vect_arch_end_aux #the initial vector of the following iteration must be the final 
      #vector of the previous iteration to compare it with the final vector returned by the swap2 in the next
      #iteration.
      vect_arch_end_aux <- vect_arch_end2 #final vector returned by the swap2 function.
    }
  }
  
  #Actual alpha coefficients for the optimal vector of archetypoids:
  zs <- x_gvv[,vect_arch_end_aux] 
  zs <- as.matrix(zs)
  alphas_def <- matrix(0, nrow = numArchoid, ncol = n)
  for (j in 1:n) {
    alphas_def[, j] = coef(nnls(zs, x_gvv[,j]))
  }
  
  resid_def <- zs[1:(nrow(zs) - 1),] %*% alphas_def - x_gvv[1:(nrow(x_gvv) - 1),]
  
  return(list(cases = vect_arch_end_aux, rss = rss_aux, archet_ini = vect_arch_ini, 
              alphas = alphas_def, resid = resid_def))
}

#' Further swap function of archetypoid analysis when looking one single archetypoid
#' @param vect_arch_end Current vector of archetypoids.
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss Current residual sum of squares of the vector of archetypoids. 
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @return Optimal vector of archetypoids.
#' @noRd
swap2_k1_norm_frob <- function(vect_arch_end, vect_arch_ini, rss, huge=200, numArchoid, x_gvv, n){
  
  vect_arch_ini_aux <- vect_arch_ini
  vect_arch_end_aux <- vect_arch_end
  
  rss_aux <- rss
  vect_arch_end2 <- c()
  
  while (vect_arch_end_aux != vect_arch_ini_aux) {
    
    rss1 <- c()
    
    setpossibles <- setdiff(1:n,vect_arch_end_aux)
    
    for (i in setpossibles) {
      zs <- x_gvv[,i] 
      zs <- as.matrix(zs)
      alphas <- matrix(0, nrow = numArchoid, ncol = n)
      for (j in 1:n) {
        alphas[, j] = coef(nnls(zs, x_gvv[,j]))
      }
      
      resid <- zs[1:(nrow(zs) - 1),] %*% alphas - x_gvv[1:(nrow(x_gvv) - 1),]
      rss1[i] <- frobenius_norm(resid) / n
      #resid <- zs %*% alphas - x_gvv
      #rss1[i] <- max(svd(resid)$d) / n
      
      if (rss1[i] < rss_aux) {
        rss_aux <- rss1[i]
        vect_arch_end2 = i
      }
    }
    
    if (is.null(vect_arch_end2)) { 
      vect_arch_end_aux <- vect_arch_end_aux
      vect_arch_ini_aux <- vect_arch_end_aux 
    }else{
      vect_arch_ini_aux <- vect_arch_end_aux 
      vect_arch_end_aux <- vect_arch_end2 
    }
  }
  
  #Actual alpha coefficients for the optimal vector of archetypoids:
  zs <- x_gvv[,vect_arch_end_aux] 
  zs <- as.matrix(zs)
  alphas_def <- matrix(0, nrow = numArchoid, ncol = n)
  for (j in 1:n) {
    alphas_def[, j] = coef(nnls(zs, x_gvv[,j]))
  }
  
  resid_def <- zs[1:(nrow(zs) - 1),] %*% alphas_def - x_gvv[1:(nrow(x_gvv) - 1),]
  
  return(list(cases = vect_arch_end_aux, rss = rss_aux, archet_ini = vect_arch_ini, 
              alphas = alphas_def, resid = resid_def))
}