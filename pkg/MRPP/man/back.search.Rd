\name{mrppBVS}
\Rdversion{1.1}
\alias{mrppBVS}
\alias{mrppBVS.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backward variable selection for the MRPP
}
\description{
\code{mrppBVS} performs backward variable selection for the MRPP. \code{mrppBVS.test} nests the backward variable selection within an outer layer of permutations to assess its significance. 
}
\usage{
mrppBVS(y,permutedTrt, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, size.in=0L, stepwise=FALSE, niter=Inf, verbose=TRUE, ...)
mrppBVS.test(y,permutedTrt, Bperm=nperms.permutedTrt(permutedTrt), 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, alpha.del=0, 
	 size.in=1L, stepwise=FALSE, verbose=TRUE, niter=Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	Numeric data matrix, with columns being variables and rows being observations.
}
  \item{permutedTrt}{
	list of matrices specifying random treatment assignments, computed from \code{\link{permuteTrt}}. 
}
  \item{importance}{
	Either \code{'dp.dw'} or \code{'p.dd.dw'}. The importance measure. \code{'dp.dw'} means the derivative of approximate p-value with respect to weights. \code{'p.dd.dw'} means the permutation p-value using the derivative of Euclidean distance with respect to weights as the new distance measure. 
}
  \item{alpha.in}{
	Desired minimally allowable level of importance for selected variables. By default, this is set to 0 if \code{importance='dp.dw'}; and 0.1 if \code{importance='p.dd.dw'}.
}
  \item{alpha.del}{
	Desired minimally allowable MRPP p-value for excluded variables. If this is set to zero, then even if the MRPP p-value for the deleted variables becomes small, the iterations continue until any other exit criteria is met.
}
  \item{size.in}{
	Minimally allowable number of selected variables. If the number of variables to be kept in the next iteration falls below this number, the iteration terminates.
}
  \item{stepwise}{
	Logical. If \code{TRUE}, then only the variables with least importance are removed. Otherwise, all variables with importance less than \code{alpha.in} are removed.
}
  \item{niter}{
	Maximum number of iterations. 
}
  \item{verbose}{
	Logical or numeric scalar. Print messages every \code{verbose} iteration. \code{TRUE} is treated as 1. If \code{FALSE}, no messages are printed.
}
%  \item{cperm.mat}{
%	Similar to \code{permutedTrt}, except that this is the permutation indices matrix for the second treatment group. 
%}
  \item{Bperm}{
	Number of outer permutations. 
}
  \item{\dots}{
	Additional arguments passed to \code{mrpp.test.dist}, \code{get.dp.dw.kde} or \code{get.p.dd.dw}.
}
}
\details{
\code{mrppBVS} performs backward variable selection for the MRPP. It starts from all \code{ncol(y)} response variables, and compute the variable importance measure specified by \code{importance}. In each iteration, it deletes a set of least important variables. When \code{stepwise=TRUE}, the variables with importance measure equal to the largest variable importance measures will delted. When \code{stepwise=FALSE}, variables with importance measure not smaller than \code{alpha.in} will be deleted. This process is iterated until the all remaining variables have importance measure not larger than \code{alpha.in}, \emph{or} the MRPP p-value for the deleted variables is not larger than \code{alpha.del}, \emph{or} the number of remaining variables will be smaller than \code{size.in} in the next iteration, \emph{or} the number of remaining variable is zero. 

\code{mrppBVS.test} performs a permutation test with the test statistic being the MRPP p-value of selected variable from the last iteration of the \code{mrppBVS} result. For simplicity, the current implementation assumes that outer permutations (that performs the final test) will be a subset of the inner permutations (that performs variable selection) specified by \code{permutedTrt}. Therefore, \code{Bperm} will be no larger than \code{ncol(permutedTrt[[1L]])}.
}
\value{
For \code{mrppBVS}, the result is a list with each component being list of :
  \item{iter }{The iteration number (starting from 0);}
  \item{var.idx }{The indices of remaining variables at this iteration, in the order of variable importance;}
  \item{influence }{The variable importance measures (multiplied by the number of permutations, if \code{importance} is \code{dp.dw}) for the remaining variables;}
  \item{p.value }{The raw MRPP p-value for the remaining variables;}
  \item{deleted.p.value }{The raw MRPP p-value for the excluded variables.}
In addition, the result also has the attribute \code{'parameter'}, which is a list of \code{importance, alpha.in, alpha.del, size.in, stepwise, niter} of input parameters, and the attribute \code{'status'}, which is a character vector giving the exit status(es).


\code{mrppBVS.test} returns an \code{htest} object with the last element of the result from \code{mrppBVS}, and 
	\item{statistic }{The observed test statistic, i.e., the raw MRPP p-value from the selectd variables.}
	\item{all.statistics }{All test statistics, i.e., the raw MRPP p-values from the selected variables for each permutation. }
	\item{p.value }{The final permutation p-value.}
	\item{parameter }{The same as the \code{parameter} attribute of \code{mrppBVS} result, converted to numeric.}
	\item{data.name }{The character name of the data.}
	\item{method }{A string of test method.}
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse
Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review)
}
\author{
Long Qu
}
\note{
%%%  ~~further notes~~
In the results, the variable importance measures are multiplied by the number of permutations, if \code{importance} is \code{dp.dw}. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{mrpp.test}}, \code{\link{get.p.dd.dw}}, \code{\link{get.dp.dw.kde}}, \code{\link{permuteTrt}}
}
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
mrpp.test(dist(x), trt, permutedTrt=pmat, wtmethod=0 ) ## MRPP for all variables
mrppBVS(x, pmat) ## backward selection of variables

\dontrun{
mrppBVS.test(x, pmat) ## nested permutation test with backward variable selection
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{htest}
\keyword{ nonparametric }
\keyword{ multivariate }
