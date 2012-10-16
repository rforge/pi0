\name{smrpp.defaultSpar}
\Rdversion{1.1}
\alias{smrpp.defaultSpar}
\title{
Default Smoothing parameters for Sparse Weighted MRPP Test 
}
\description{
Compute a sequence of default smoothing parameters covering the minimum (sparsest) and a user specified maximum, plus infinity (densest).
}
\usage{
smrpp.defaultSpar(dp.dw, max.ratio=2, nspar=100L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp.dw}{A vector or matrix of derivative of MRPP p-value to the weights. If it is a matrix, each row is treated as a permutation, as returned from \code{\link{get.dp.dw.kde}} with \code{test=TRUE}.
}
  \item{max.ratio}{A positive numeric scalar. Let \code{k} be the ratio of maximum weight to the minimum weight.
 When spar increases, \code{k}  decreases. 
 Whenever \code{k} drops below \code{max.ratio}, the corresponding \code{spar} will be treated as \code{Inf} (resulting in \code{k=1})
}
  \item{nspar}{An integer no smaller than 3, requesting the number of smoothing parameters to be generated.
}
}

\value{
\code{smrpp.defaultSpar} returns a vector or matrix of smoothing parameters. 

When \code{dp.dw} is a vector, the result is also a vector. 
When \code{dp.dw} is a matrix, a separate sequence of smoothing parameters will be generated for each row. 

In each row, the smoothing parameters are sorted. The smallest (i.e., the first) one is given by 
\deqn{\frac{\imath^{(2)}-\imath^{(1)}}{2R/m_1},}{(See reference for the formula)}
where \eqn{\imath^{(1)}}{i^(1)} and \eqn{\imath^{(2)}}{i^(2)} are the smallest and the 2nd smallest \emph{unique} \code{dp.dw} values for the row, 
\eqn{R}{R} is the total number of variables, and \eqn{m_1}{m_1} is the number of \code{dp.dw} values that tie with the smallest value in the row. 
For each row of the result, the last element is always \code{Inf}. The second last element is determined by the \code{max.ratio} as described above. 
Except for the last one, the spacing between successive smoothing parameters are equal on the log10 scale. 

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
 \code{\link{smrpp.test}},  \code{\link{smrpp.penWt}}, 
}
\examples{
set.seed(2340L)
dp.dw=runif(10L, -1, 1)  # some fake derivatives
smrpp.defaultSpar(dp.dw, max.ratio=2, nspar=5L) # a vector of smoothing parameters

set.seed(2340L)
dp.dw=matrix(runif(100L, -1, 1),10L,byrow=TRUE)  # some fake derivatives
smrpp.defaultSpar(dp.dw, max.ratio=2, nspar=5L) # a vector of smoothing parameters for each row of dp.dw

smrpp.defaultSpar(dp.dw, max.ratio=2, nspar=3L)[,1L] # getting the minimum of smoothing parameters 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smoothing}
