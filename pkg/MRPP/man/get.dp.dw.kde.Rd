\name{MRPP importance}
\Rdversion{1.1}
\alias{get.dp.dw.kde}
\alias{get.p.dd.dw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation-based variable importance measure
}
\description{
\code{get.dp.dw.kde} computes variable importance measure of partial derivative of continuous approximation to MRPP p-values with respect to hypothetical weights. 
\code{get.p.dd.dw} computes variable importance measure by using partial derivative of distance function with respect to hypothetical weights as the new distance function in MRPP
}
\usage{
get.dp.dw.kde(y, permutedTrt, r = seq_len(ncol(y)), test = FALSE, distObj = dist(y), 
	mrpp.stats = mrpp.test.dist(distObj, permutedTrt = permutedTrt,
	wtmethod=wtmethod[1])$all.statistics, bw = bw.mse.pdf.asym(mrpp.stats), 
	wtmethod=0, scale=1, standardized=FALSE)
get.p.dd.dw(y, permutedTrt, r = seq_len(ncol(y)), wtmethod=0, eps=1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The response data matrix with each row being a data point and each column being a variable. 
}
  \item{permutedTrt}{A list matricies that list all random treatment assignments to be used in the MRPP test. 
}
  \item{r}{A positive integer vector, specifying variables for which variable importance measure is requested. 
}
  \item{test}{\code{TRUE} or \code{FALSE}. If \code{FALSE}, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. If \code{TRUE}, each permutation in \code{permutedTrt} will be treated as if it is the original treatment assignment, useful for doing permutation tests using the computed variable importance meausre as the test statistic.
}
  \item{distObj}{The numeric vector computed from \code{dist(y)}, i.e., the lower triangle of the Euclidean distance matrix. 
}
  \item{mrpp.stats}{A vector of all observed MRPP test statistics. 
}
  \item{bw}{A numeric scalar giving the bandwidth used for kernel density estimation of \code{mrpp.stats}. This might be positive infinity. 
}
%  \item{cperm.mat}{}
  \item{wtmethod}{0 or 1. If 0, the treatment group weight for MRPP statistic will be sample size minus 1. If 1, the sample size of the treatment group.
}
  \item{scale}{A numeric scalar. The result will be multiplied by this scale factor. 
}
  \item{standardized}{In case that \code{bw} is positive infinity, \code{standardized} indicates if the results will be divided by the standard deviation of the permuted statistics.
}
%  \item{min.wts}{}
  \item{eps}{A small positive tolerance, the same as described in \code{\link{mrpp.test}}. 
}
}
\details{
\code{get.dp.dw.kde} first uses kernel density estimator to approximate the distribution of \code{mrpp.stats} and obtain a continous approximate p-value for the MRPP test. Then it computes the partial derivative of the approximate p-value with respect to the hypothetical weight of each variable, evaluated at weight 1. 
That is, if the partial derivative is negatively large, then the corresponding variable will be more important, because increasing its weight will quickly decreasing the approximate MRPP p-value. Vice versa. 

\code{get.p.dd.dw} first computes the partial derivative of the Euclidian distance measure with respect to the hypothetical weights for each variable. This is then used as the distance function for the variable to be used in MRPP. The resulting MRPP p-value is returned as the variable importance. The smaller, the more important.
}
\value{
For \code{get.dp.dw.kde}, when \code{test} is \code{FALSE}, the result is an \code{length(r)}-vector of importance measure (multiplied by \code{scale}). 
When \code{test} is \code{TRUE}, the result is an \code{B}-by-\code{length(r)} matrix of importance measure, where \code{B} is the number of permutations in \code{permutedTrt}. The (\code{i,j})th element is the importance measure for the \code{r[j]}th variable when the \code{i}th column of each component of \code{permutedTrt} is treated as the original treatment assignment. 

For \code{get.p.dd.dw}, the result is an \code{length(r)}-vector of importance measure. 
}
\note{
When \code{test} is \code{TRUE} and when the number of permutations is large, \code{get.dp.dw.kde} will consume a large amount of memory. 
}

\references{
Long Qu, Dan Nettleton, Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review). 
}
\author{
Long Qu
}

\seealso{\code{\link{mrpp.test}}, \code{\link{permuteTrt}}, \code{\link{bw.mse.pdf.asym}}, \code{\link{mrppBVS.test}}, \code{\link{mrppBVS}}
}
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
ir1=get.dp.dw.kde(x, pmat)
ir2=get.p.dd.dw(x, pmat)
cor(ir1, ir2, method='s')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}
\keyword{nonparametric}