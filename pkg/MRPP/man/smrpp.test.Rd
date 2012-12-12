\name{smrpp.test}
\Rdversion{1.1}
\alias{smrpp.test}
\alias{smrpp.test.default}
\alias{smrpp.test.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse Weighted MRPP Test for One-way Design
}
\description{
K-sample sparse MRPP test with penalization. Variable selection is through penalized weighting. The choice of tuning parameter is through minimizing the MRPP p-value, with weighted Euclidean distance. Such a minimum p-value is then used as the final test statistic. The overal significance is through an outer layer of permutations. 
}
\usage{
smrpp.test(y, \dots)
\method{smrpp.test}{default}(y, trt, B=nparts(table(trt)), permutedTrt, wtmethod=0, 
	outerStat=c('WDISCO 1/F','WMRPP P'), eps=1e-8, spar, verbose=TRUE, \dots)
\method{smrpp.test}{formula}(y, data, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An N*R data matrix, with rows representing samples and columns representing variables.
}
  \item{trt}{A vector of treatment assignments. If missing, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. As least one of \code{trt} or \code{permutedTrt} need to be given. 
}
  \item{B}{A positive integer number of permutations requested. This will only be used when \code{permutedTrt} is missing. This will be overwritten by the \code{ncol(permutedTrt)}.
}
  \item{permutedTrt}{An optional permutation index matrix. If missing, this will be computed as \code{permuteTrt(trt,B)}. See \code{\link{permuteTrt}}. 
}
  \item{wtmethod}{0 or 1, where 0 stands for weighting each treatment group by sample size - 1, and 1 stands for weighting by sample size. 
}
  \item{outerStat}{A character value specifying which test statistic to use for the outer permutation. 'WDISCO 1/F' uses the weighted DISCO inverse F-statistic. 'WMRPP P' uses weighted MRPP raw p-value. 
}
  \item{eps}{A small non-negative number, differences below which among permuted test statistics are treated as equal.
}
  \item{spar}{A positive numeric vector of smoothing parameters, over which the minimum p-value will be searched for. An \code{Inf} will be padded to \code{spar}, if it is not already there.
}
  \item{data}{A data frame in which variables in the formula can be found.}
  \item{verbose}{
	Logical or numeric scalar. Print messages every \code{verbose} iteration. \code{TRUE} is treated as 1. If \code{FALSE}, no messages are printed.
}
  \item{\dots}{Additional arguments passed to the methods. 
}
}


\value{

An \code{htest} object with the following components:
\item{statistic }{The observed test statistic (the minimum MRPP p-value).}
\item{all.statistics }{All permuted test statistics.}
\item{weights }{The selected weights for each dimension of the response variable.}
\item{p.value }{The permutation p-value.}
\item{parameter }{A vector of the number of permutations, the weighting method and the chosen smoothing parameter.}
\item{data.name }{The character name of the data.}
%\item{.Random.seed }{ The random number seed for reproducibility.} 
\item{method }{A string of test method.}
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
 \code{\link{smrppInitSpar}}, 
 \code{\link{smrpp.penWt}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
dat=data.frame(x, trt)
fmla=as.formula(sprintf('cbind(\%s)~trt',paste('X',1:5,sep='',collapse=',')))

#\dontrun{
smrpp.test(x, trt, permutedTrt=pmat, wtmethod=0 ) 
smrpp.test(fmla, dat, permutedTrt=pmat, wtmethod=0 )
#}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}