\encoding{latin1}
\name{summary.groc}
\alias{summary.groc}
\alias{print.groc}
\title{Summary and Print Methods for groc objects}
\description{
  Summary and print methods for \code{groc} objects.
}
\usage{
\method{summary}{groc}(object, what = "validation",
        digits = 4, print.gap = 2, \dots)
\method{print}{groc}(x, \dots)
}
\arguments{
  \item{x, object}{a \code{groc} object}
  \item{what}{character, only \code{"validation"} for the moment}
  \item{digits}{integer.  Minimum number of significant digits in the
    output.  Default is 4.}
  \item{print.gap}{Integer.  Gap between coloumns of the printed tables.}
  \item{\dots}{Other arguments sent to underlying methods.}
}
\details{
  If \code{what} is \code{"validation"}, the cross-validated PRESS,
  RPEMAD and RMSEPs (if
  available) are given.
}
\value{
  \code{print.groc} return the object invisibly.
}
\references{
Martin Bilodeau, Pierre Lafaye de Micheaux, Smail Mahdi (2015), The R
Package groc for Generalized Regression on Orthogonal Components,
\emph{Journal of Statistical Software}, 65(1), 1-29, \cr \url{http://www.jstatsoft.org/v65/i01/}
}
\author{P. Lafaye de Micheaux}
\seealso{\code{\link{groc}}, \code{\link{grocCrossval}}}
\examples{
data("yarn",package="pls")
yarn.groc <- groc(density ~ NIR, 6, data = yarn)
yarn.cv <- grocCrossval(yarn.groc, segments = 10)
print(yarn.groc)
summary(yarn.cv)
}
\keyword{regression}
\keyword{multivariate}
