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

\dontrun{
########################
# Codes for Example 1  #
########################
require("groc")
data("wood")
out <- groc(y ~ x1 + x2 + x3 + x4 + x5, ncomp = 1, data = wood, 
             D = corrob, method = "lts")
corrob(wood$y, fitted(out)) ^ 2
plot(out)

########################
# Codes for Example 2  #
########################
data("trees")
out <- groc(Volume ~ Height + Girth, ncomp = 1, D = spearman, 
             method = "s", data = trees)
cor(trees$Volume, fitted(out)) ^ 2
plot(out$T, trees$Volume, xlab = "First latent variable",
     ylab = "Volume", pch = 20)
lines(sort(out$T), fitted(out)[order(out$T)])
out <- boxcox(Volume ~ Height + Girth, data = trees, 
              lambda = seq(-0.5, 0.5, length = 100), plotit = FALSE)
lambda <- out$x[which.max(out$y)]
out <- lm(Volume ^ lambda ~ Height + Girth, data = trees)
cor(trees$Volume, fitted(out)^(1/lambda)) ^ 2

########################
# Codes for Example 3  #
########################
data("wood")
plsr.out <- plsr(y ~ x1 + x2 + x3 + x4 + x5, data = wood)
groc.out <- groc(y ~ x1 + x2 + x3 + x4 + x5, data = wood)
apply(abs((fitted(plsr.out) - fitted(groc.out)) / 
          fitted(plsr.out)), 3, max) * 100

########################
# Codes for Example 4  #
########################
set.seed(1)
n <- 200
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
y <- x1 * x2 + rnorm(n, 0, sqrt(.04))
data <- data.frame(x1 = x1, x2 = x2, y = y)
plsr.out <- plsr(y ~ x1 + x2, data = data)
groc.out <- groc(y ~ x1 + x2, D = dcov, method = "s", data = data)
plsr.v <- crossval(plsr.out, segment.type = "consecutive")
groc.v <- grocCrossval(groc.out, segment.type = "consecutive")
groc.v$validation$PRESS
plsr.v$validation$PRESS
gam.data <- data.frame(y = y, t1 = groc.out$T[, 1], t2 = groc.out$T[, 2])
gam.out <- gam(y ~ s(t1) + s(t2), data = gam.data)
par(mfrow = c(1, 2))
plot(gam.out)
par(mfrow = c(1, 1))
PRESS <- 0
for(i in 1 : 10){
  data.in <- data[-(((i - 1) * 20 + 1) : (i * 20)), ]
  data.out <- data[((i - 1) * 20 + 1) : (i * 20), ]
  ppr.out <- ppr(y ~ x1 + x2, nterms = 2, optlevel = 3, data = data.in)
  PRESS <- PRESS + sum((predict(ppr.out, newdata = data.out)-data.out$y) ^ 2)
}
PRESS

########################
# Codes for Example 5  #
########################
data("yarn")
dim(yarn$NIR)
n <- nrow(yarn)
system.time(plsr.out <- plsr(density ~ NIR, ncomp = n - 2, data = yarn))
system.time(groc.out <- groc(density ~ NIR, Nc = 20, ncomp = n - 2, data = yarn))
max(abs((fitted(plsr.out) - fitted(groc.out)) / fitted(plsr.out))) * 100
plsr.v <- crossval(plsr.out, segments = n, trace = FALSE)
plsr.v$validation$PRESS
groc.v <- grocCrossval(groc.out, segments = n, trace = FALSE)
groc.v$validation$PRESS
groc.v$validation$PREMAD

########################
# Codes for Example 6  #
########################
data("prim7")
prim7.out <- groc(X1 ~ ., ncomp = 3, D = dcov, method = "s", data = prim7)
prim7.out$R
pca <- princomp(~ ., data = as.data.frame(prim7[, -1]))
prim7.pca <- data.frame(X1 = prim7$X1, scores = pca$scores)
prim7.pca.out <- groc(X1 ~ ., ncomp = 3, D = dcov, method = "s", 
                       data = prim7.pca)
pca$loadings %*% prim7.pca.out$R
groc.v <- grocCrossval(prim7.out, segment.type = "consecutive")
groc.v$validation$PRESS
plsr.out <- plsr(X1 ~ ., ncomp = 3, data = prim7)
plsr.v <- crossval(plsr.out, segment.type = "consecutive")
plsr.v$validation$PRESS
PRESS <- 0
for(i in 1 : 10){
  data.in <- prim7[-(((i - 1) * 50 + 1) : (i * 50)), ]
  data.out <- prim7[((i - 1) * 50 + 1) : (i * 50), ]
  ppr.out <- ppr(X1 ~ ., nterms = 3, optlevel = 3, data = data.in)
  PRESS <- PRESS + sum((predict(ppr.out, newdata = data.out) - data.out$X1) ^ 2)
}
PRESS

########################
# Codes for Example 7  #
########################
n <- 50 ; B <- 30
mat.cor <- matrix(0, nrow = B, ncol = 3) ; mat.time <- matrix(0, nrow = B, ncol = 3)
for (i in 1:B) {
 X <- matrix(runif(n * 5, -1, 1), ncol = 5)
 A <- matrix(runif(n * 50, -1, 1), nrow = 5)
 y <- (X[,1] + X[,2])^2 + (X[,1] + 5 * X[,2])^2 + rnorm(n)
 X <- cbind(X, X %*% A)
 D <- data.frame(X = X, y = y)
 mat.time[i,1] <- system.time(out1 <- plsr(y ~ X, , ncomp = 2, data = D))[1]
 mat.time[i,2] <- system.time(out2 <- ppr(y ~ X, , nterms = 2, data = D))[1]
 mat.time[i,3] <- system.time(out3 <- groc(y ~ X, D = dcov, method = "s", ncomp = 2, data = D))[1]
 mat.cor[i,] <- cor(y, cbind(fitted(out1)[,,2], fitted(out2), fitted(out3)[,,2]))
}
colMeans(mat.cor)
colMeans(mat.time)

########################
# Codes for Example 8  #
########################
data("oliveoil")
n <- nrow(oliveoil)
plsr.out <- plsr(sensory ~ chemical, data = oliveoil, method = "simpls")
groc.out <- groc(sensory ~ chemical, data = oliveoil)
max(abs((fitted(plsr.out) - fitted(groc.out)) / fitted(plsr.out))) * 100
groc.v <- grocCrossval(groc.out, segments = n)
groc.v$validation$PRESS
colMeans(groc.v$validation$PRESS)
Y <- oliveoil$sensory
for (j in 1 : ncol(Y)) print(cor(Y[, j], fitted(groc.out)[, j, 2]))

########################
# Codes for Example 9  #
########################
require("ppls")
data("cookie")
X <- as.matrix(log(cookie[1 : 40, 51 : 651]))
Y <- as.matrix(cookie[1 : 40, 701 : 704])
X <- X[, 2 : 601] - X[, 1 : 600]
data <- data.frame(Y = I(Y), X = I(X))
n <- nrow(data)
q <- ncol(Y)
xl <- "Wavelength index"
yl <- "First differences of log(1/reflectance)"
matplot(1:ncol(X), t(X), lty = 1, xlab = xl, ylab = yl, type = "l")
out1 <- plsr(Y ~ X, ncomp = n - 2, data = data)
cv <- crossval(out1, segments = n)
cv.mean <- colMeans(cv$validation$PRESS)
plot(cv.mean, xlab = "h", ylab = "Average PRESS", pch = 20)
h <- 3
for (j in 1 : q) print(cor(Y[, j], fitted(out1)[, j, h]))
set.seed(1)
out2 <- groc(Y ~ X, ncomp = h, data = data, plsrob = TRUE)
for (j in 1 : q) print(corrob(Y[, j], fitted(out2)[, j, h]))
plot(out2)

########################
# Codes for Example 10 #
########################
set.seed(2)
n <- 30
t1 <- sort(runif(n, -1, 1))
y <- t1 + rnorm(n, mean = 0, sd = .05)
y[c(14, 15, 16)] <- y[c(14, 15, 16)] + .5
data <- data.frame(x1 = t1, x2 = 2 * t1, x3 = -1.5 * t1, y = y)
out <- groc(y ~ x1 + x2 + x3, ncomp = 1, data = data, plsrob = TRUE)
tau <- scaleTau2(residuals(out), mu.too = TRUE)
std.res <- scale(residuals(out), center = tau[1], scale = tau[2])
index <- which(abs(std.res)>3)
prm.res <- read.table("prmresid.txt")
plot(t1, y, pch = 20)
matlines(t1, cbind(t1,fitted(out), y - prm.res), lty = 1 : 3)
legend(.4, -.5 , legend = c("true model","groc", "prm"), lty = 1 : 3)
text(t1[index], y[index], index, cex = .8, pos = 3)

########################
# Codes for Example 11 #
########################
data("pulpfiber")
X <- as.matrix(pulpfiber[, 1:4])
Y <- as.matrix(pulpfiber[, 5:8])
data <- data.frame(X = I(X), Y = I(Y))
set.seed(55481)
out.rob <- groc(Y ~ X, data = data, plsrob = TRUE)
plot(out.rob, cex = .6)
out.simpls <- groc(Y ~ X, data = data)
cv.rob <- grocCrossval(out.rob,segment.type = "consecutive")
PREMAD.rob <- cv.rob$validation$PREMAD[,4]
PREMAD.rob
cv.simpls <- grocCrossval(out.simpls,segment.type = "consecutive")
PREMAD.simpls <- cv.simpls$validation$PREMAD[,4]
PREMAD.simpls
(PREMAD.rob - PREMAD.simpls) / PREMAD.simpls * 100
}

}
\keyword{distribution} % Probability Distributions and Random Numbers
\keyword{htest} % Statistical Inference
