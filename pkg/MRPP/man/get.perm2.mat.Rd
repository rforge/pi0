\name{get.perm2.mat}
\Rdversion{1.1}
\alias{get.perm2.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute permutation index matrix
}
\description{
Compute permutation index matrix for 2-sample design
}
\usage{
get.perm2.mat(trt, B = 100L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{A vector or factor with exactly two levels, indicating the treatment assignment.

}
  \item{B}{A positive integer indicating the number of permutations required
}
}
\details{
This is a wrapper to \code{\link[pi0]{combn2R}}.
}
\value{
An integer matrix, such that the number of rows is the \code{min(table(trt))}, i.e., the number of samples in the smaller treatment group, and  the number of columns is the number of permutations produced. Note that, when \code{B} is larger than the number of possible treeatment assignments, the number of columns of the resulting matrix could be less than \code{B}.

The result also has the \code{.Random.seed} attribute, which is the random seed at the start of the call. 
The other attribute \code{sample.method} is the same as in \code{\link[pi0]{combn2R}}.

In addition, the first column of the result will always be the (unpermuted) original treatment assignment.
}
\author{
Long Qu
}
\seealso{
\code{\link[pi0]{combn2R}}
}
\examples{
set.seed(2340)
get.perm2.mat(gl(2,3))[,]
get.perm2.mat(gl(2,3),10L)[,]
str(get.perm2.mat(gl(2,5)))
}
\keyword{utilities}
\keyword{iteration}
