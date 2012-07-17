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
\method{mrpp.test}{dist}(y, trt, B = nparts(attr(y, "Size"), table(trt)), perm.mat, wtmethod=0, eps=1e-8, \dots)
\method{mrpp.test}{formula}(y, data, \dots)
\method{mrpp.test}{default}(y, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A N*P data matrix, \code{\link{dist}} object of the data matrix, or a \code{\link{formula}}. 
}
  \item{trt}{A vector of treatment assignments. 
}
  \item{B}{A positive integer number of permutations requested. This would be overwritten by \code{ncol(perm.mat)} if \code{perm.mat} is non-missing.
}
  \item{perm.mat}{An optional permutation index matrix. See \code{\link{get.perm.mat}}. 
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
\item{.Random.seed }{ The random number seed for reproducibility.} 
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
 \code{\link{get.perm.mat}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts(length(trt), table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(1,seed=1032940L) # init seed
pmat=get.perm.mat(trt, 5e2L)		## use 500 random permutations
\dontrun{
pmat=get.perm.mat(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
mrpp.test(dist(x), trt, perm.mat=pmat, wtmethod=0 ) 
}
mrpp.test(dist(x), trt, perm.mat=pmat, wtmethod=0 ) 
mrpp.test(x, trt, perm.mat=pmat, wtmethod=0 ) 
dat=data.frame(x, trt)
fmla=as.formula(sprintf('cbind(\%s)~trt',paste0('X',1:5,collapse=',')))
mrpp.test(fmla, dat, perm.mat=pmat, wtmethod=0 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
