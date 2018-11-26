#' archetypoids_funct_multiv_robust_internal
#' 
#' @rdname archetypoids_funct_multiv_robust_internal
#'
#' @description 
#' These three functions are called by the archetypoid algorithm 
#' (\code{\link{archetypoids_funct_multiv_robust}}). All these functions are therefore not 
#' solely used.
#' 
#' @author 
#' Irene Epifanio
#' 
#' @seealso 
#' \code{\link{archetypoids_funct_multiv_robust}}
#' 
#' @references 
#' Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
#' with application to financial time series analysis, 2018, submitted,
#' \url{https://arxiv.org/abs/1810.00919}
#'                  
#' General swap function of the functional multivariate robust archetypoid analysis 
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss_arch_ini Residual sum of squares associated to \code{vect_arch_ini}.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
#' @return Optimal vector of archetypoids.
#' @noRd
swap_funct_multiv_robust <- function(vect_arch_ini, rss_arch_ini, huge = 200, numArchoid, 
                                  x_gvv, n, PM, prob, nbasis, nvars){
  
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
      rss1[i] <- frobenius_norm_funct_multiv_robust(resid, PM, prob, nbasis, nvars) / n
      
      if (rss1[i] < rss) { 
        rss <- rss1[i]
        vect_arch_end = c(i,vect_arch_ini[-l])
      }
    }
  }
  
  if (numArchoid == 1) {
    result <- swap2_k1_funct_multiv_robust(vect_arch_end, vect_arch_ini, rss, huge, 
                                        numArchoid, x_gvv, n, PM, prob, nbasis, nvars)
  }else{
    result <- swap2_funct_multiv_robust(vect_arch_end, vect_arch_ini, rss, huge, 
                                     numArchoid, x_gvv, n, PM, prob, nbasis, nvars)
  } 
  
  return(result)
}

#' Further swap function of the functional multivariate robust archetypoid analysis 
#' @param vect_arch_end Current vector of archetypoids.
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss Current residual sum of squares of the vector of archetypoids.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
#' @return Optimal vector of archetypoids.
#' @noRd
swap2_funct_multiv_robust <- function(vect_arch_end, vect_arch_ini, rss, huge = 200, numArchoid, 
                                   x_gvv, n, PM, prob, nbasis, nvars){
  
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
        rss1[i] <- frobenius_norm_funct_multiv_robust(resid, PM, prob, nbasis, nvars) / n
        
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

#' Further swap function of the functional multivariate robust archetypoid analysis 
#' when looking one single archetypoid
#' @param vect_arch_end Current vector of archetypoids.
#' @param vect_arch_ini Initial vector of archetypoids.
#' @param rss Current residual sum of squares of the vector of archetypoids.
#' @param huge Penalization added to solve the convex least squares problems. 
#' @param numArchoid Number of archetypoids to compute.
#' @param x_gvv Original data matrix with the added penalization.
#' @param n Number of observations in data.
#' @param PM Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.
#' @param prob Probability with values in [0,1].
#' @param nbasis Number of basis.
#' @param nvars Number of variables.
#' @return Optimal vector of archetypoids.
#' @noRd
swap2_k1_funct_multiv_robust <- function(vect_arch_end, vect_arch_ini, rss, huge=200, 
                                      numArchoid, x_gvv, n, PM, prob, nbasis, nvars){
  
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
      rss1[i] <- frobenius_norm_funct_multiv_robust(resid, PM, prob, nbasis, nvars) / n
      
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