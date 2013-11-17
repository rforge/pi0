\name{MRPP-package}
\Rdversion{1.1}
\alias{MRPP-package}
\alias{MRPP}
\docType{package}
\title{
Multiresponse permutation procedure and its variable importance and selection. 
}
\description{
Implements multiresponse permutation procedure (MRPP), variable importance measure for MRPP, and backward variable selection procedures for the MRPP, applicable to high dimensional data. 
}
\details{
\tabular{ll}{
Package: \tab MRPP\cr
Type: \tab Package\cr
Version: \tab 1.0-321\cr
Date: \tab 2013-11-17\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
\itemize{
\item{\code{\link{mrpp.test}}} {The MRPP, a permutation test applicable to high dimensional data.}
\item{\code{\link{get.dp.dw.kde}}} {Variable importance measure: approximate derivative of p-values to weights.}
\item{\code{\link{get.p.dd.dw}}} {Variable importance measure: permutation p-value using the derivative of Euclidean distance with respect to weights as the distance measure.}
\item{\code{\link{mrppBVS}}} {Backward variable selection for the MRPP.}
\item{\code{\link{mrppBVS.test}}} {Backward variable selection nested within outer permutations.}
\item{\code{\link{smrpp.test}}} {Sparse Weighted MRPP within outer permutations.}
\item{\code{\link{smrpp.penWt}}} {Sparse Weights given MRPP importance measures from \code{\link{get.dp.dw.kde}}.}
\item{\code{\link{permuteTrt}}} {Generate random treatment assignments.}
\item{\code{\link{FR2permvec}}} {Conversions among decimal factoradic number, radix factoradic number, and permutation vector.}
\item{\code{\link{HSEL.bigz}}} {Generate a random sample of size m from a \emph{large} population of size n.}
\item{\code{\link{bw.mse.pdf.asym}}} {Choose a bandwidth for kernel density estimation by minimizing asymptotic mean squared error of density estimate at the first data point.}
}
}
\author{
Long Qu

Maintainer: Long Qu <long.qu@wright.edu>
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse
Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review)
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
mrpp.test(dist(x), trt, permutedTrt=pmat, wtmethod=0 ) ## MRPP for all variables

\dontrun{
mrppBVS(x, pmat) ## backward selection of variables
mrppBVS.test(x, pmat) ## nested permutation test with backward variable selection
smrpp.test(x, trt, permutedTrt=pmat, wtmethod=0 ) 
}

}
