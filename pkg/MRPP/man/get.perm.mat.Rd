\name{permuteTrt}
\Rdversion{1.1}
\alias{permuteTrt}
\alias{nparts}
\alias{nperms.permutedTrt}
\alias{ntrt.permutedTrt}
\alias{trt.permutedTrt}
\title{
Compute permutations of treatment assignments in a one-way design. 
}
\description{
Enumerate or sample permutations of treatment assignments in a one-way design. 
}
\usage{
permuteTrt(trt, B = 100, idxOnly = FALSE)
nparts(n)
nperms.permutedTrt(permutedTrt)
ntrt.permutedTrt(permutedTrt)
trt.permutedTrt(permutedTrt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{A vector of treatment labels
}
  \item{B}{A positive integer of the number of permutations requested. 
}
  \item{idxOnly}{A logical scalar, indicating whether only factoradic indices are computed. This has no effect when \code{B} is not smaller than \code{nparts(table(trt))}. See the Value section for details.}
  \item{n}{The sample size in each treatment group.}
  \item{permutedTrt}{A \code{'permutedTrt'} object from \code{permutedTrt} call.}
}
\details{
The \code{nparts} functions computes the total number of partitions, with sizes being specified by \code{n}. 

Denote the totality by \code{N}, i.e., \code{length(trt)} or \code{sum(n)}. When \code{B} is no smaller than \code{nparts(table(trt))}, all partitions are returned by \code{permuteTrt} and \code{B} is set to \code{nparts(table(trt))}. 
Otherwise, a sample \code{B} random permutations will be selection from a total of \code{N!} possible permutations (without replacement). 
}
\value{
For \code{permuteTrt}, it returns S3 object of class \code{'permutedTrt'}, which is a list \code{ANS} of \code{length(unique(trt))} treatment assignment matrices. When \code{idxOnly} is \code{FALSE}, each matrix has \code{B} columns, where \code{B} may be smaller than requested, such that \code{ANS[[i]][,b]} contains sorted individual indices that are allocated to treatment \code{i} under the \code{b}th random permutation. Note that, here \code{b} should not be interpreted as the factoradic number. That is, the column locations of the matrices have no special meaning. However, an exception is that the first column (i.e., the first permutation) is always equivalent to the original treatment assignment. In other words, it is guaranteed that \code{setequal( split(seq(length(trt)),trt), lapply(ANS, '[', , 1L))} should always be \code{TRUE}. In addition, the \code{'idx'} attribute of the result is set ot \code{NA_character_}. 

When \code{idxOnly} is \code{TRUE}, each matrix has only one column, corresponding to the original treatment assignment. The permutations being sampled represented by the correspoding decimal factoradic numbers (stored as character) returned in the \code{'idx'} attribute. In this case, the length of the \code{'idx'} attribute is \code{B} and the first element is always \code{"0"}, corresponding to no permutations.

\code{nparts} returns the total number of non-equivalent treatment assignments. 
\code{nperms.permutedTrt} returns the total number of permutations. 
\code{ntrt.permutedTrt} returns a named vector of sample sizes. 
\code{trt.permutedTrt} recovers the \code{trt} vector passed to \code{permuteTrt}.
}
\note{
When \code{B} is less than \code{nparts(table(trt))}, ideally one could randomly sample \code{B} partitions from a total of \code{nparts(table(trt))} without replacement, subject to the first partition being equivalent to the original treatment assignment. For the sake of simplicity, this is \emph{not} implemented here. The random sample of size \code{B} is generated from a total of \code{N!} permutations without replacement. Therefore, it is possible for some of the permutations sampled corresponding to equivalent partitions. However, for the purpose of permutation tests, the type I error rate is still under control. 
}
\author{
Long Qu
}
\seealso{
\code{\link{combn}}, \code{\link[gmp]{urand.bigz}}, \code{\link{HSEL.bigz}}, \code{\link{sample1}}, \code{\link{dec2permvec}}
}
\examples{
set.seed(2340)
trt=gl(2,10)
nparts(table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 1e3L)		## use 1000 random permutations
\dontrun{
pmat=permuteTrt(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{iteration}
\keyword{arith}