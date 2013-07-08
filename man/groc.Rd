\name{groc}
\encoding{utf8}
\alias{groc}
\alias{groc.default}
\title{groc method}
\description{
Generalized regression on orthogonal components.
}
\usage{
\method{groc}{default}(formula, ncomp, data, subset, na.action, plsrob =
                 FALSE, method = c("lm", "lo", "s", "lts"), D = NULL,
                 gamma = 0.75, Nc = 10, Ng = 20, scale = FALSE, Cpp =
                 TRUE, model = TRUE, x = FALSE, y = FALSE, ...)
groc(\dots)
}
\arguments{
 \item{formula}{a model formula.  Most of the \code{lm} formula
    constructs are supported.  See below.}
  \item{ncomp}{the number of components (orthogonal components) to include in the model.}
  \item{data}{an optional data frame with the data to fit the model from.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain missing values.}
  \item{plsrob}{logical. If \code{TRUE}, we use the \code{D=covrob} measure of
    dependence with the least trimmed squares method="lts".}
  \item{method}{character giving the name of the method to use. The
    user can supply his own function. The methods available are linear
    models, "lm", local polynomials, "lo", smoothing splines, "s", and least trimmed squares, "lts".}
  \item{D}{function with two arguments, each one being a vector, which
    measures the dependence between two variables using n observations from them. If \code{NULL}, the covariance measure will be used. The user can supply his own function.}
  \item{gamma}{parameter used with the option \code{plsrob=TRUE}. It defines the quantile used to compute the "lts" regression. The default \code{gamma=0.75} gives a breakdown of 25\% for a good compromise between robustness and efficiency. The value \code{gamma=0.5} gives the maximal breakdown of 50\%.}
  \item{Nc}{Integer, Number of cycles in the grid algorithm.}
  \item{Ng}{Integer, Number of points for the grid in the grid algorithm.}
  \item{scale}{Logical, Should we scale the data.}
  \item{Cpp}{Logical, if \code{TRUE} this function will use a C++
    implementation of the grid algorithm. The \code{FALSE} value should not be
    used, unless to get a better understanding of the grid algorithm
    or to compare the speed of computation between R and C++ versions of
  this algorithm}
  \item{model}{a logical.  If \code{TRUE}, the model frame is returned.}
  \item{x}{a logical.  If \code{TRUE}, the model matrix is returned.}
  \item{y}{a logical.  If \code{TRUE}, the response is returned.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
%\details{
%TODO
%}
\value{
  \item{Y}{vector or matrix of responses.}
  \item{fitted.values}{an array of fitted values.}
  \item{residuals}{residuals}
  \item{T}{a matrix of orthogonal components (scores). Each column corresponds to a component.}
  \item{R}{a matrix of directions (loadings). Each column is a direction used to obtain the corresponding component (scores).}
  \item{Gobjects}{contain the objects produced by the fit of the responses on the orthogonal components.}
  \item{Hobjects}{contain the objects produced by the "lts" fit of each deflated predictors on the orthogonal components. \code{Hobjects} are produced when \code{plsrob=TRUE}.}
  \item{B}{matrix of coefficients produced by the "lm" fit of  each deflated predictors on the last component. \code{B} is produced when \code{plsrob=FALSE}.}
  \item{Xmeans}{a vector of means of the X variables.}
  \item{Ymeans}{a vector of means of the Y variables.}
  \item{D}{Dependence measure used.}
  \item{V}{a matrix whose columns contain the right singular vectors of
    the data. Computed in the preprocessing to principal component scores when the number of
observations is less than the number of predictors.}
  \item{dnnames}{dimnames of 'fitted.values'}
  \item{ncomp}{the number of components used in the modelling.}
  \item{method}{the method used.}
  \item{scale}{Logical. \code{TRUE} if the responses have been scaled.}
  \item{call}{the function call.}
  \item{terms}{the model terms.}
  \item{plsrob}{Logical. If \code{plsrob=TRUE}, a robust partial least squares fit.}
  \item{model}{if \code{model=TRUE}, the model frame.}

}
%\references{
%}
\author{Martin Bilodeau (\email{bilodeau@dms.umontreal.ca}) and Pierre Lafaye de Micheaux (\email{lafaye@dms.umontreal.ca}) and
  Smail Mahdi (\email{smail.mahdi@cavehill.uwi.edu})
}
\examples{
data(wood)
out <- groc(y ~ x1+x2+x3+x4+x5, ncomp=1, data=wood,
D=corrob, method="lts")
corrob(wood$y,fitted(out))^2
plot(out)

}
\keyword{distribution} % Probability Distributions and Random Numbers
\keyword{htest} % Statistical Inference
