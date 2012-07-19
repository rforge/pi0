\name{permuteTrt}
\Rdversion{1.1}
\alias{permuteTrt}
\alias{nparts}
\title{
Compute permutations of treatment assignments in a one-way design. 
}
\description{
Compute permutations of treatment assignments in a one-way design. 
}
\usage{
permuteTrt(trt, B = 100)
nparts(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{A vector of treatment labels
}
  \item{B}{A positive integer of the number of permutations requested. 
}
%  \item{N}{The total number of observations.}
  \item{n}{The sample size in each treatment group.}
}
\details{
The \code{nparts} functions computes the total number of partitions, with sizes being specified by \code{n}. 

Denote the totality by \code{N}, i.e., \code{length(trt)} or \code{sum(n)}. When \code{B} is no smaller than \code{nparts(table(trt))}, all partitions are returned by \code{permuteTrt} and \code{B} is set to \code{nparts(table(trt))}. 
Otherwise, a sample \code{B} random permutations will be selection from a total of \code{N!} possible permutations (without replacement). 
}
\value{
For \code{permuteTrt}, a list \code{ANS} of \code{length(unique(trt))} treatment assignment matrices, with each matrix having \code{B} columns, where \code{B} may be smaller than requested, such that ANS[[i]][,b] contains sorted individual indices that are allocated to treatment \code{i} under the \code{b}th random permutation. Note that, here \code{b} should not be interpreted as the factoradic number. That is, the column locations of the matrices have no special meaning. However, an exception is that the first column (i.e., the first permutation) is always equivalent to the original treatment assignment. In other words, it is guaranteed that \code{setequal( split(seq(length(trt)),trt), lapply(ANS, '[', , 1L))} should always be \code{TRUE}. 

for \code{nparts}, the total number of non-equivalent treatment assignments. 
}
\note{
When \code{B} is less than \code{nparts(table(trt))}, ideally one could randomly sample \code{B} partitions from a total of \code{nparts(table(trt))} without replacement, subject to the first partition being equivalent to the original treatment assignment. For the sake of simplicity, this is \emph{not} implemented here. The random sample of size \code{B} is generated from a total of \code{N!} permutations without replacement. Therefore, it is possible for some of the permutations sampled corresponding to equivalent partitions. However, for the purpose of permutation tests, the type I error rate is still under control. 
}
\author{
Long Qu
}
\seealso{
\code{\link{combn}}, \code{\link[gmp]{urand.bigz}}, \code{\link{HSEL.bigz}}, \code{\link{sample1}}
}
\examples{
set.seed(2340)
trt=gl(2,10)
nparts(table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(1,seed=1032940L) # init seed
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