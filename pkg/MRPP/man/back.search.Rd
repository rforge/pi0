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
back.search(y, permutedTrt, verbose = TRUE, niter = Inf, 
	    importance = c("dp.dw", "p.dd.dw"), 
	    alpha.in, alpha.del = 0, stepwise = FALSE, 
	    size.in=0L, ...)
bsmrpp(y,permutedTrt, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, alpha.del=0, stepwise=FALSE, 
	 size.in=1L, Bperm=ncol(permutedTrt), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	Numeric data matrix, with columns being variables and rows being observations.
}
  \item{permutedTrt}{
	A permutation of indices of the first treatment group. Each column represent a permutation.
}
  \item{verbose}{
	Logical or numeric scalar. Print messages every \code{verbose} iteration. \code{TRUE} is treated as 1. If \code{FALSE}, no messages are printed.
}
  \item{niter}{
	Maximum number of iterations. 
}
  \item{importance}{
	Either \code{'dp.dw'} or \code{'p.dd.dw'}. The importance measure. \code{'dp.dw'} means the derivative of approximate p-value with respect to weights. \code{'p.dd.dw'} means the permutation p-value using the derivative of Euclidean distance with respect to weights as the new distance measure. 
}
  \item{alpha.in}{
	Desired minimally allowable level of importance for selected variables. By default, this is set to 0 if \code{importance='dp.dw'}; and 0.1 if \code{importance='p.dd.dw'}.
}
  \item{alpha.del}{
	Desired minimally allowable MRPP p-value for excluded variables. 
}
  \item{stepwise}{
	Logical. If \code{TRUE}, then only the variables with least importance are removed. Otherwise, all variables with importance less than \code{alpha.in} are removed.
}
  \item{size.in}{
	Minimally allowable number of selected variables. 
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{back.search} returns a list, each element corresponding to results in each iteration. Each element is a list by itself, with elements
\item{iter}{Iteration index;}
\item{var.idx}{Variable index, in the order of variable importance;}
\item{influence}{Sorted ariable importance;}
\item{p.value}{MRPP p-value for selected variables;}
\item{deleted.p.value}{MRPP p-value for excluded variables.}

\code{bsmrpp} returns the last element of the result from \code{back.search}, except that
\item{raw.p.value}{replaces the original p-value from \code{back.search};}
\item{adj.p.value}{the actual p-value (adjusted) computed by nesting variable selection in permutations;}
\item{permuted.p.values}{a vector of raw.p.values, one for each outer permutation.}
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. RELATIVE VARIABLE IMPORTANCE AND BACKWARD VARIABLE SELECTION FOR THE MULTIRESPONSE
PERMUTATION PROCEDURE, WITH APPLICATIONS TO HIGH DIMENSIONAL GENOMIC DATA. (submitted)
}
\author{
Long Qu
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mrpp.test}}, \code{\link{get.p.dd.dw}}, \code{\link{get.dp.dw.kde}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
