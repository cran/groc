\name{corrob}
\alias{corrob}
\title{Robust correlation measure}
\description{
Compute robust estimates of the correlation between two variables using the Orthogonalized Gnanadesikan-Kettenring pairwise estimator. 
}
\usage{
corrob(t, u)
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
This function uses the \code{\link{covRob}} function from the \pkg{robust} package.
}
\value{
Value of the robust correlation.
}
\references{
  Jiahui Wang, Ruben Zamar, Alfio Marazzi, Victor Yohai, Matias
  Salibian-Barrera, Ricardo Maronna, Eric Zivot, David Rocke, Doug
  Martin, Martin Maechler and Kjell Konis. (2013). robust: Robust
  Library. R package version 0.4-11.
  \url{https://CRAN.R-project.org/package=robust}
}
\author{
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@unsw.edu.au})
}
\seealso{
\code{\link{covrob}}, \code{\link{dcov}}
}
\examples{
data(stackloss)
corrob(stackloss$Air.Flow,stackloss$Water.Temp)
}
\keyword{robust}
