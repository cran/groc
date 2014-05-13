\name{predict.groc}
\alias{predict.groc}
\title{Predict Method for groc
}
\description{
Prediction for groc models.  New responses or scores are
  predicted using a fitted model and a new matrix of observations.
}
\usage{
\method{predict}{groc}(object, newdata, ncomp = object$ncomp, na.action = na.pass, ...)
}
\arguments{
  \item{object}{
a \code{groc} object.  The fitted model
}
  \item{newdata}{
a data frame.  The new data.  If missing, the training
    data is used.
}
  \item{ncomp}{
vector of positive integers.  The components to use in
    the prediction.}
  \item{na.action}{
function determining what should be done with missing
    values in \code{newdata}. By default, nothing is done.
}
  \item{\dots}{
further arguments.  Currently not used
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A three dimensional array of
  predicted response values is returned.  The dimensions correspond to
  the observations, the response variables and the model sizes,
  respectively.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
%%  ~~who you are~~
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@dms.umontreal.ca})
}

\seealso{
\code{\link{plot.groc}}
}
\examples{
data(wood)
out <- groc(y ~ x1+x2+x3+x4+x5, ncomp=1, data=wood,D=corrob, method="lts")
predict(out)

newdata<- data.frame(x1= 0.5, x2=0.1, x3=0.4, x4=0.5, x5=0.8)
predict(out,newdata)
}
\keyword{regression}
\keyword{multivariate}
