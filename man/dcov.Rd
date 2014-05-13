\name{dcov}
\alias{dcov}
\title{Distance covariance matrix.
}
\description{
Compute the distance covariance measure of Szekely, Rizzo, and Bakirov
(2007) between two samples.
}
\usage{
dcov(x, y, Cpp = TRUE)
}
\arguments{
  \item{x}{
data of first sample
}
  \item{y}{
data of second sample
}
  \item{Cpp}{
logical. If TRUE (the default), computations are performed using a C version of the code.
  }
}
\details{
See \pkg{energy}.
}
\value{
  returns the sample distance covariance.
}
\references{
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), 
 Measuring and Testing Dependence by Correlation of Distances, 
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \url{http://dx.doi.org/10.1214/009053607000000505}
}
\author{
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@dms.umontreal.ca})
}
\seealso{
\code{\link{covrob}}, \code{\link{corrob}}
}
\examples{
data(stackloss)
dcov(stackloss$Air.Flow,stackloss$Water.Temp)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance covariance }
