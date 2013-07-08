\name{covrob}
\alias{covrob}
\title{Robust covariance measure}
\description{
Compute robust estimates of the covariance between two variables using
the robust tau estimate of univariate scale, as proposed by Maronna and Zamar (2002).
}
\usage{
covrob(t, u)
}
\arguments{
  \item{t}{
a numeric vector containing the data for the fisrt variable.
}
  \item{u}{
a numeric vector containing the data for the second variable.
}
}
\details{
This function uses the \code{\link{scaleTau2}} function from the \pkg{robustbase} package.
}
\value{
Value of the robust covariance.
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.
}
\author{
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@dms.umontreal.ca})
}
\seealso{
\code{\link{corrob}}, \code{\link{dcov}}
}
\examples{
data(stackloss)
covrob(stackloss$Air.Flow,stackloss$Water.Temp)
}
\keyword{robust}
