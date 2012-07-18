\name{back.search}
\Rdversion{1.1}
\alias{back.search}
\alias{bsmrpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backward variable selection for the MRPP
}
\description{
\code{back.search} performs backward variable selection for the MRPP. \code{bsmrpp} nests the backward variable selection within an outer layer of permutations to assess its significance. 
}
\usage{
back.search(y,permutedTrt, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, size.in=0L, stepwise=FALSE, niter=Inf, verbose=TRUE, ...)
bsmrpp(y,permutedTrt, Bperm=ncol(permutedTrt[[1L]]), 
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
\code{back.search} performs backward variable selection for the MRPP. It starts from all \code{ncol(y)} response variables, and compute the variable importance measure specified by \code{importance}. In each iteration, it deletes a set of least important variables. When \code{stepwise=TRUE}, the variables with importance measure equal to the largest variable importance measures will delted. When \code{stepwise=FALSE}, variables with importance measure not smaller than \code{alpha.in} will be deleted. This process is iterated until the all remaining variables have importance measure not larger than \code{alpha.in}, \emph{or} the MRPP p-value for the deleted variables is not larger than \code{alpha.del}, \emph{or} the number of remaining variables will be smaller than \code{size.in} in the next iteration, \emph{or} the number of remaining variable is zero. 

%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{back.search}, the result is a list with each component being list of :
  \item{iter }{The iteration number (starting from 0);}
  \item{var.idx }{The indices of remaining variables at this iteration, in the order of variable importance;}
  \item{influence }{The variable importance measures for the remaining variables;}
  \item{p.value }{The raw MRPP p-value for the remaining variables;}
  \item{deleted.p.value }{The raw MRPP p-value for the excluded variables.}
In addition, the result also has attribute \code{'setting'}, which is a list of \code{importance, alpha.in, alpha.del, size.in, stepwise, niter} of input parameters. 


\code{bsmrpp} returns the last element of the result from \code{back.search}, except that
  \item{raw.p.value }{replaces the original p-value from \code{back.search};}
  \item{adj.p.value }{the actual p-value (adjusted) computed by nesting variable selection in permutations;}
  \item{permuted.p.values }{a vector of raw.p.values, one for each outer permutation.}
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse
Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review)
}
\author{
Long Qu
}
%\note{
%%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{mrpp.test}}, \code{\link{get.p.dd.dw}}, \code{\link{get.dp.dw.kde}}, \code{\link{permuteTrt}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
