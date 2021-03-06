% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/int_prod_mat_sq.R
\name{int_prod_mat_sq}
\alias{int_prod_mat_sq}
\title{Squared interior product between matrices}
\usage{
int_prod_mat_sq(m)
}
\arguments{
\item{m}{Data matrix.}
}
\value{
Data matrix.
}
\description{
Helper function to compute the robust Frobenius norm.
}
\examples{
mat <- matrix(1:4, nrow = 2)
int_prod_mat_sq(mat)

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
