\name{groc.fit}
\alias{groc.fit}
\title{Fitting a groc model}
\description{
  Fits a groc model with the grid algorithm.
}
\usage{
groc.fit(X, Y, ncomp = min(nrow(X) - 1, ncol(X)), D = NULL, gamma =
                 0.75, method = NULL, plsrob = FALSE, Nc = 10, Ng = 20,
                 scale = FALSE, Cpp = TRUE, stripped = FALSE, maxiter =
                 100, sp = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a matrix of predictors.  \code{NA}s and \code{Inf}s are not allowed.
}
  \item{Y}{
a vector or matrix of responses.  \code{NA}s and \code{Inf}s are not allowed.
}
  \item{ncomp}{
the number of components to be used in the modelling.
}
  \item{D}{
Dependence measure.}
  \item{gamma}{
    Used to set the breakdown value when \code{method="lts"}.
  }
  \item{method}{
the method to be used. Currently only 'lm', 'lo', 's', and 'lts'.
}
  \item{plsrob}{
Logical. If \code{TRUE}, the function sets \code{D=covrov} and \code{method="lts"} for a robust partial least squares fit.
}
\item{Nc}{
Integer. Number of cycles in the grid algorithm}
  \item{Ng}{
Integer. Number of points for the grid in the grid algorithm.}
  \item{scale}{
Logical. If \code{TRUE} the responses are scaled.
}
  \item{Cpp}{
Logical. If \code{TRUE}, compuations are performed in a faster way using
a C code.}
  \item{stripped}{
logical.  If \code{TRUE} the calculations are stripped
    as much as possible for speed; this is meant for use with
    cross-validation or simulations when only the coefficients are
    needed.  Defaults to \code{FALSE}.
}
  \item{maxiter}{
Integer. Maximal number of iterations in the grid algorithm. Used only when there are more than one response.
}
\item{sp}{
      A vector of smoothing parameters can be provided here.
          Smoothing parameters must be supplied in the order that the
          smooth terms appear in the model formula. Negative elements
          indicate that the parameter should be estimated, and hence a
          mixture of fixed and estimated parameters is
	  possible. 'length(sp)' should be equal to 'ncomp' and
	  corresponds 
	  to the number of underlying smoothing parameters.
  }
  \item{\dots}{
other arguments.  Currently ignored.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{Y}{data used as response.}
  \item{fitted.values}{an array of fitted values. Its element [i,j,k] is the fitted value for observation i, response j, and when  k components are used.}
  \item{residuals}{an array of regression residuals.  It has the same
    dimensions as \code{fitted.values}.}
  \item{T}{a matrix of orthogonal components (scores). Each column corresponds to a component.}
  \item{R}{a matrix of directions (loadings). Each column is a direction used to obtain the corresponding component (scores).}
  \item{Gobjects}{contain the objects produced by the fit of the responses on the orthogonal components.}
  \item{Hobjects}{contain the objects produced by the "lts" fit of each deflated predictors on the orthogonal components. \code{Hobjects} are produced when \code{plsrob=TRUE}.}
  \item{B}{matrix of coefficients produced by the "lm" fit of  each deflated predictors on the last component. \code{B} is produced when \code{plsrob=FALSE}.}
  \item{Xmeans}{a vector of means of the X variables.}
  \item{Ymeans}{a vector of means of the Y variables.}
  \item{D}{Dependence measure used.}
  \item{V}{a matrix whose columns contain the right singular vectors of the data.
Computed in the preprocessing to principal component scores when the number of
observations is less than the number of predictors.}
  \item{dnnames}{dimnames of 'fitted.values'}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de
Micheaux (\email{lafaye@dms.umontreal.ca})}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{multivariate}% __ONLY ONE__ keyword per line
