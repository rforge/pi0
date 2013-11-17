\name{mrpp.test.dist}
\Rdversion{1.1}
\alias{mrpp.test.dist}
\alias{mrpp.test.default}
\alias{mrpp.test.formula}
\alias{mrpp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MRPP test for one-way design
}
\description{
K-sample MRPP test
}
\usage{
mrpp.test(y, \dots)
\method{mrpp.test}{dist}(y, trt, B = as.integer(min(nparts(table(trt)), 1e4L)), permutedTrt, wtmethod=0, eps=1e-8, \dots)
\method{mrpp.test}{formula}(y, data, \dots)
\method{mrpp.test}{default}(y, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A N*P data matrix, \code{\link{dist}} object of the data matrix, or a \code{\link{formula}}. 
}
  \item{trt}{A vector of treatment assignments. If missing, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. As least one of \code{trt} or \code{permutedTrt} need to be given. 
}
  \item{B}{A positive integer number of permutations requested. This will only be used when \code{permutedTrt} is missing. This will be overwritten by the \code{ncol(permutedTrt)}.
}
  \item{permutedTrt}{An optional permutation index matrix. If missing, this will be computed as \code{permuteTrt(trt,B)}. See \code{\link{permuteTrt}}. 
}
  \item{wtmethod}{0 or 1, where 0 stands for weighting each treatment group by sample size - 1, and 1 stands for weighting by sample size. 
}
  \item{eps}{A small non-negative number, differences below which among permuted test statistics are treated as equal.
}
  \item{data}{A data frame in which variables in the formula can be found.}
  \item{\dots}{Additional arguments passed to the methods. 
}
}

\value{
An \code{htest} object with the following components:
\item{statistic }{The observed test statistic.}
\item{all.statistics }{All permuted test statistics.}
\item{p.value }{The permutation p-value.}
\item{parameter }{A vector of the number of permutations and the weighting method.}
\item{data.name }{The character name of the data.}
%\item{.Random.seed }{ The random number seed for reproducibility.} 
\item{method }{A string of test method.}
}
\references{
Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
}
\author{
Long Qu
}
%%\note{
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{permuteTrt}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
\dontrun{
pmat=permuteTrt(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
}
### Distance matrix interface
mrpp.test(dist(x), trt, permutedTrt=pmat, wtmethod=0 ) 

### Data matrix interface
mrpp.test(x, trt, permutedTrt=pmat, wtmethod=0 ) 

### Formula interface
dat=data.frame(x, trt)
fmla=as.formula(sprintf('cbind(\%s)~trt',paste('X',1:5,sep='',collapse=',')))
mrpp.test(fmla, dat, permutedTrt=pmat, wtmethod=0 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}