\name{get.perm.mat}
\Rdversion{1.1}
\alias{get.perm.mat}
\alias{nparts}
\title{
Compute permutation matrix
}
\description{
Compute permutation matrix based on design
}
\usage{
get.perm.mat(trt, B = 100)
nparts(N, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{A vector of treatment labels
}
  \item{B}{A positive integer of the number of permutations requested. 
}
  \item{N}{The total number of observations.}
  \item{n}{The sample size in each treatment group.}
}
\details{
When \code{B} is larger than \code{nparts(N, table(trt))}, all partitions are returned and \code{B} is adjusted. 
Otherwise, a sample \code{B} random permutations will be selection from a total of \code{N!}. 
}
\value{
For \code{get.perm.mat}, a list of \code{n}-by-\code{B} permutation index matrices; 
for \code{nparts}, the total number of non-equivalent treatment assignments. 
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
nparts(length(trt), table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(1,seed=1032940L) # init seed
pmat=get.perm.mat(trt, 1e3L)		## use 1000 random permutations
\dontrun{
pmat=get.perm.mat(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{iteration}
