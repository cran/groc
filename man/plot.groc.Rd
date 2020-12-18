\name{plot.groc}
\alias{plot.groc}
\title{Plot groc objects.
}
\description{
A function to plot groc objects.
}
\usage{
\method{plot}{groc}(x, h=x$ncomp, cex=0.8, ...)
}
\arguments{
  \item{x}{
A groc object.
}
  \item{h}{
Number of components in the model.
  }
  \item{cex}{
Character expansion factor for point labels.}
  \item{\dots}{
Further arguments passed to internal \code{plot} function.}
}
\details{
If \code{plsrob=FALSE}, a plot of robust Mahalanobis distances for residuals versus robust Mahalanobis distances for components. Useful for identification of good points, vertical outliers, good and bad leverage points.

If \code{plsrob=TRUE}, the previous plot is done with another similar plot of classical Mahalanobis distances to compare the identification
of the various type of points obtained by classical or robust partial least squares. 
}
\references{
Martin Bilodeau, Pierre Lafaye de Micheaux, Smail Mahdi (2015), The R
Package groc for Generalized Regression on Orthogonal Components,
\emph{Journal of Statistical Software}, 65(1), 1-29, \cr \url{https://www.jstatsoft.org/v65/i01/}
}
\author{
%%  ~~who you are~~
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@unsw.edu.au})
}


\examples{
## This example takes some time:
\dontrun{
data("pulpfiber",package="robustbase")
X <- as.matrix(pulpfiber[, 1:4])
Y <- as.matrix(pulpfiber[, 5:8])
data <- data.frame(X=I(X), Y=I(Y))
set.seed(55481)
out.rob <- groc(Y ~ X, data=data, plsrob=TRUE)
plot(out.rob, cex=.6)
}
}
\keyword{regresion}
\keyword{multivariate}
