\name{smrpp.penWt}
\Rdversion{1.1}
\alias{smrpp.penWt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing Penalized Weights for Sparse Weighted MRPP Test for One-way Design
}
\description{
This function solves the penalized optimization problem by minimizing an approximation of weighted MRPP p-values subject to constraints on the size of weight adjustment, total weight and nonnegativity of weights.
}
\usage{
smrpp.penWt(dp.dw, spar, simplify=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp.dw}{A vector or matrix of derivative of MRPP p-value to the weights. If it is a matrix, each row is treated as a permutation, as returned from \code{\link{get.dp.dw.kde}} with \code{test=TRUE}.
}
  \item{spar}{A nonnegative numeric vector of smoothing parameters, or a matrix. If it is a matrix, \code{dp.dw} needs also be a matrix with the same number of rows. All computation will be performed for each row separately.  
}
  \item{simplify}{A logical scalar. If \code{TRUE}, the penalized weight will be passed to \code{drop()} before being returned.
}
}
\details{
Let \code{lambda} be a smoothing parameter under consideration (e.g., one value of \code{spar}) and \code{Ir} be a vector of derivatives of MRPP p-value to the weights (e.g., a row of \code{dp.dw}). 
This function finds \code{w} that minimizes \code{ sum(Ir * w) + lambda * sum(w * w) } subject to the constraints \code{all(w>=0)} and \code{sum(w)==R}, where \code{R} is the total number of variables.
}

\value{
When \code{dp.dw} is a vector, it is converted to a single-row matrix. The result is then a 3D array with \code{dim} being \code{c(L, B, R)}, 
where \code{L} is the number of smoothing parameters, \code{B} is the number of rows of \code{dp.dw}, and \code{R} is the number of variables. 
If \code{simplify=TRUE}, this array is passed to \code{drop()} before being returned. 
}

\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. Relative Variable Importance and Variable Selection for the Multiresponse Permutation Procedure, with Applications to High Dimensional Genomic Data.
}
\author{
Long Qu
}
%%\note{
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{smrpp.defaultSpar}},
 \code{\link{smrpp.test}}
}
\examples{
set.seed(2340L)
dp.dw=runif(10L, -1, 1)  # some fake derivatives
spar=smrpp.defaultSpar(dp.dw, max.ratio=2, nspar=5L) # a vector of smoothing parameters
smrpp.penWt(dp.dw, spar, TRUE) # a matrix of weights, with rows of decreasing sparseness
smrpp.penWt(dp.dw, spar, FALSE) # the same result without dropping the 2nd array dimension
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{htest}
%\keyword{multivariate}
\keyword{smoothing}
\keyword{optimize}