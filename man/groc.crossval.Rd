\name{grocCrossval}
\alias{grocCrossval}
\title{Cross-validation of groc models
}
\description{
  A \dQuote{stand alone} cross-validation function for \code{groc} objects.
}
\usage{
grocCrossval(object, segments = 10, segment.type = c("random",
      "consecutive","interleaved"), length.seg, trace = 15, ...)
}
\arguments{
  \item{object}{
a \code{groc} object; the regression to cross-validate.
}
  \item{segments}{
    the number of segments to use, or a list with segments (see below).
  }
  \item{segment.type}{
the type of segments to use.
}
  \item{length.seg}{
Positive integer.  The length of the segments to use.
}
  \item{trace}{
if \code{TRUE}, tracing is turned on.  If numeric, it
    denotes a time limit (in seconds).  If the estimated total time of
    the cross-validation exceeds this limit, tracing is turned on.
}
  \item{\dots}{
additional arguments, sent to the underlying fit function.
}
}
\details{
  This function performs cross-validation on a model fit by \code{groc}.
  It can handle models such as \code{groc(Y ~ X, \dots)}.

  Note that to use \code{grocCrossval}, the data \emph{must} be specified
  with a \code{data} argument when fitting \code{object}.

  If \code{segments} is a list, the arguments \code{segment.type} and
  \code{length.seg} are ignored.  The elements of the list should be
  integer vectors specifying the indices of the segments.  

  Otherwise, segments of type \code{segment.type} are generated.  How
  many segments to generate is selected by specifying the number of
  segments in \code{segments}, or giving the segment length in
  \code{length.seg}.  If both are specified, \code{segments} is
  ignored.

  When tracing is turned on, the segment number is printed for each segment.
}
\value{
  The supplied \code{object} is returned, with an additional component
  \code{validation}, which is a list with components
  \item{method}{equals \code{"CV"} for cross-validation.}
  \item{pred}{an array with the cross-validated predictions.}
  \item{PRESS}{a matrix of PRESS values for models with 1, \ldots,
    \code{ncomp} components.  Each row corresponds to one response variable.}
  \item{PREMAD}{a matrix of PREMAD values for models with 1, \ldots,
    \code{ncomp} components.  Each row corresponds to one response variable.}
  \item{RMSEP}{a matrix of sqrt(PRESS/nobj) values for models with 1, \ldots,
    \code{ncomp} components.  Each row corresponds to one response variable.}
  \item{segments}{the list of segments used in the cross-validation.}
  \item{ncomp}{the number of components.}
}
\references{
Martin Bilodeau, Pierre Lafaye de Micheaux, Smail Mahdi (2015), The R
Package groc for Generalized Regression on Orthogonal Components,
\emph{Journal of Statistical Software}, 65(1), 1-29, \cr \url{https://www.jstatsoft.org/v65/i01/}
}
\author{
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@unsw.edu.au})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(yarn,package="pls")
yarn.groc <- groc(density ~ NIR, 6, data = yarn)
yarn.cv <- grocCrossval(yarn.groc, segments = 10)

yarn.cv$validation$PRESS
yarn.cv$validation$PREMAD
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{multivariate}% __ONLY ONE__ keyword per line
