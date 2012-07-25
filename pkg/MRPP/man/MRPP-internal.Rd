\name{internal}
\Rdversion{1.1}
\alias{get.p.L1dEuc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions not be used by the end user
}
\description{
\code{get.p.L1dEuc} computes variable importance measure by using (dimension specific L1 distance divided by Euclidean distance across all dimensions) as the new distance function in MRPP
}
\usage{
get.p.L1dEuc(y, permutedTrt, r = seq_len(ncol(y)), wtmethod=0, eps=1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The response data matrix with each row being a data point and each column being a variable. 
}
  \item{permutedTrt}{A list matricies that list all random treatment assignments to be used in the MRPP test. 
}
  \item{r}{A positive integer vector, specifying variables for which variable importance measure is requested. 
}
  \item{wtmethod}{0 or 1. If 0, the treatment group weight for MRPP statistic will be sample size minus 1. If 1, the sample size of the treatment group.
}
  \item{eps}{A small positive tolerance, the same as described in \code{\link{mrpp.test}}. 
}
}
\details{
\code{get.p.L1dEuc} is very similar to \code{get.p.dd.dw}. It first computes the L1 distance for each variable separately, then it divides this by the overall Euclidean distance across all dimensions. This is then used as the distance function for the variable to be used in MRPP. The resulting MRPP p-value is returned as the variable importance. The smaller, the more important.
}
\value{
For \code{get.p.L1dEuc}, the result is an \code{length(r)}-vector of importance measure. 
}

\references{
Long Qu, Dan Nettleton, Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review). 
}
\author{
Long Qu
}

\seealso{\code{\link{get.p.dd.dw}}
}
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
ir3=MRPP:::get.p.L1dEuc(x, pmat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{internal}