\encoding{latin1}
\name{model.frame.groc}
\alias{model.frame.groc}
\alias{model.matrix.groc}
\title{Extract Information From a Fitted groc Model}
\description{
  Functions to extract information from \code{groc} objects:  the model frame, the model
  matrix.
}
\usage{
\method{model.matrix}{groc}(object, \dots)
\method{model.frame}{groc}(formula, \dots)
}
\arguments{
  \item{object, formula}{a \code{groc} object.  The fitted model.}
  \item{\dots}{other arguments sent to underlying functions.}
}
\details{

  \code{model.frame.groc} returns the model frame; i.e. a data frame with
  all variables necessary to generate the model matrix.  See
  \code{\link[stats]{model.frame}} for details.
  
\code{model.matrix.groc} returns the (possibly coded) matrix used as
  \eqn{X} in the fitting.  See \code{\link[stats]{model.matrix}} for
  details.

}
\value{
  \code{model.frame.groc} returns a data frame with
  all variables neccessary to generate the model matrix.

 \code{model.matrix.groc} returns the \eqn{X} matrix.

 }
\author{Ron Wehrens and Bjørn-Helge Mevik}
\seealso{\code{\link{coef}}, \code{\link{fitted}},
  \code{\link{residuals}}, \code{\link{model.frame}}
  }

%\examples{
%}
\keyword{regression}
\keyword{multivariate}
