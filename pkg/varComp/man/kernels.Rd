\name{kernels}
\alias{ibs}
\alias{IBS}
\alias{am}
\alias{AM}
\alias{lin0}
\alias{Lin0}
\alias{quad1}
\alias{Quad1}
\alias{Minkowski}
\alias{minkowski}
\alias{polyk}
\alias{Polyk}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Kernel functions useful for genetic associations
}
\description{
	These are the functions that might be used in computing pairwise inter-individual similarities based on their single nucleotide polymorphism (SNP) genotypes. 
}
\usage{
am(x)
AM(x)
ibs(x)
IBS(x)
lin0(x)
Lin0(x)
quad1(x)
Quad1(x)
Minkowski(x, p = 1)
minkowski(x, p = 1)
polyk(x,c=0,d=1) 
Polyk(x,c=0,d=1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A numeric matrix encoding genotypes. Each row corresponds to an individual and each column corresponds to a genetic marker. Usually, allele-counting coding is used, but others are allowed. 
}  
  \item{p}{ The exponent defining the Minkowski distance. The same as in \code{\link[stats:dist]{stats::dist}}.
}
  \item{c}{ The constant added to cross-products before raising to the power of \code{d}.
}
  \item{d}{ The exponent defining the polynomial kernel. When \code{c=0} and \code{d=1}, this is equivalent to \code{lin0}. When \code{c=1} and \code{d=2}, this is equivalent to \code{quad1}. 
}
}
\details{
These functions compute the pairwise similarities among rows of \code{x}. Lower-case versions are more useful in the formula interface to specify random genetic effects. Upper-case versions can be used to directly compute the genetic similarity matrix. 

\code{am} and \code{AM} calculate the allele-matching kernel, and \code{AM} is based on \code{\link[SPA3G:KERNEL]{SPA3G:::KERNEL}}. 

\code{ibs} and \code{IBS} compute the identity-by-descent (IBS) kernel. \code{IBS} is computed as \cr
\code{1 - as.matrix(dist(x, method='manhattan') * .5 /max(1, ncol(x)) )}. 

\code{lin0} and \code{Lin0} compute the linear kernel with zero intercept. \code{Lin0} is computed as \cr \code{normalizeTrace(tcrossprod(x)/max(1,ncol(x)))}. 

\code{quad1} and \code{Quad1} compute the quadratic kernel with offset 1. \code{Qaud1} is computed as \cr  \code{normalizeTrace((base::tcrossprod(x)+1)^2)}. 

\code{minkowski} and \code{Minkowski} compute the similarity based on the Minkowski distance. \code{Minkowski} is computed as \code{1-as.matrix(dist(x, method='minkowski', p=p)) * .5 / max(1, ncol(x))^(1/p) }. 

}
\value{
The functions starting with an upper-case letter returns an \code{n}-by-\code{n} symmetric similarity matrix, where \code{n} equals \code{nrow(x)}. The corresponding functions starting with a lower-case letter returns a matrix \code{L} such that \code{tcrossprod(L)} equals the value from their upper-case counterparts. The number of rows is \code{n}, but the number of columns is the rank of the similarity matrix. 
}
% \references{
% ~put references to the literature/web site here ~
% }
\author{
Long Qu
}
% \note{
%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cholRoot}}, \code{\link{normalizeTrace}}, \code{\link[stats:dist]{stats::dist}}, \code{\link[SPA3G:KERNEL]{SPA3G:::KERNEL}}, \code{\link{varComp}}
}
\examples{
set.seed(2345432L)
x=matrix(sample(2, 50L, replace=TRUE), 10L)  
IBS(x)
range(tcrossprod(ibs(x)) - IBS(x)  )

AM(x)
range(tcrossprod(am(x)) - AM(x)  )

Lin0(x)
range(tcrossprod(lin0(x)) - Lin0(x)  )
range(Lin0(x) - Polyk(x, 0, 1))

Quad1(x)
range(tcrossprod(quad1(x)) - Quad1(x)  )
range(Quad1(x) - Polyk(x, 1, 2))

Minkowski(x)
range(tcrossprod(minkowski(x)) - Minkowski(x)  )
range(tcrossprod(minkowski(x)) - IBS(x)  )

## Use in formulas
model.matrix(~0+ibs(x))
range(tcrossprod(model.matrix(~0+ibs(x))) - IBS(x))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% __ONLY ONE__ keyword per line
