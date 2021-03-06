% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archetypoids_funct.R
\name{archetypoids_funct}
\alias{archetypoids_funct}
\title{Archetypoid algorithm with the functional Frobenius norm}
\usage{
archetypoids_funct(numArchoid, data, huge = 200, ArchObj, PM)
}
\arguments{
\item{numArchoid}{Number of archetypoids.}

\item{data}{Data matrix. Each row corresponds to an observation and each column 
corresponds to a variable. All variables are numeric.}

\item{huge}{Penalization added to solve the convex least squares problems.}

\item{ArchObj}{The list object returned by the 
\code{\link{stepArchetypesRawData_funct}} function.}

\item{PM}{Penalty matrix obtained with \code{\link[fda]{eval.penalty}}.}
}
\value{
A list with the following elements:
\itemize{
\item cases: Final vector of archetypoids.
\item rss: Residual sum of squares corresponding to the final vector of archetypoids.
\item archet_ini: Vector of initial archetypoids.
\item alphas: Alpha coefficients for the final vector of archetypoids.
\item resid: Matrix with the residuals.
}
}
\description{
Archetypoid algorithm with the functional Frobenius norm 
to be used with functional data.
}
\examples{
\dontrun{
library(fda)
?growth
str(growth)
hgtm <- t(growth$hgtm)
# Create basis:
basis_fd <- create.bspline.basis(c(1,ncol(hgtm)), 10)
PM <- eval.penalty(basis_fd)
# Make fd object:
temp_points <- 1:ncol(hgtm)
temp_fd <- Data2fd(argvals = temp_points, y = growth$hgtm, basisobj = basis_fd)
data_archs <- t(temp_fd$coefs)

lass <- stepArchetypesRawData_funct(data = data_archs, numArch = 3, 
                                    numRep = 5, verbose = FALSE, 
                                    saveHistory = FALSE, PM)

af <- archetypoids_funct(3, data_archs, huge = 200, ArchObj = lass, PM) 
str(af)                                
}                                                          
                                                     
}
\references{
Epifanio, I., Functional archetype and archetypoid analysis, 2016. 
\emph{Computational Statistics and Data Analysis} \bold{104}, 24-34, 
\url{https://doi.org/10.1016/j.csda.2016.06.007}
}
\seealso{
\code{\link[Anthropometry]{archetypoids}}
}
\author{
Irene Epifanio
}
