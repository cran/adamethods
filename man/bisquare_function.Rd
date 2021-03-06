% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bisquare_function.R
\name{bisquare_function}
\alias{bisquare_function}
\title{Bisquare function}
\usage{
bisquare_function(resid, prob, ...)
}
\arguments{
\item{resid}{Vector of residuals, computed from the 
\eqn{m \times n} residuals data matrix.}

\item{prob}{Probability with values in [0,1].}

\item{...}{Additional possible arguments.}
}
\value{
Vector of real numbers.
}
\description{
This function belongs to the bisquare family of loss functions.
The bisquare family can better cope with extreme outliers.
}
\examples{
resid <- c(2.47, 11.85)  
bisquare_function(resid, 0.8)

}
\references{
Moliner, J. and Epifanio, I., Robust multivariate and functional archetypal analysis 
with application to financial time series analysis, 2019. 
\emph{Physica A: Statistical Mechanics and its Applications} \bold{519}, 195-208. 
\url{https://doi.org/10.1016/j.physa.2018.12.036}
}
\author{
Irene Epifanio
}
