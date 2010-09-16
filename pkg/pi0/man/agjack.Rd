\name{   agjack.pi0       }
\alias{    agjack.pi0     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Averaged generalized jackknife estimate of pi0
}
\description{
Averaged generalized jackknife estimate of pi0

}
\usage{
agjack.pi0(subtobj,mean.n=c('mean','harmean','geomean') ,
    pointpair=FALSE, trunc=TRUE,tol=1e-5)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subtobj}{      A \code{\link{subt}} object        }
  \item{mean.n}{A character: name of the function to compute the average sample sizes}
  \item{pointpair}{logical: if \code{TRUE}, then \code{\link{gjack}} is called for each pair of rows in \code{subtobj}; otherwise, \code{\link{gjack}} is called for each pair of unique average sample sizes. }
  \item{trunc}{logical, indicating if each \code{\link{gjack}} estiamte is truncated to [0,1]. }
  \item{tol}{A small tolerance number. }
}
\details{
When \code{pointpair} is \code{FALSE}, the rows in \code{subtobj} are first grouped by combination of sample sizes and the estimates are averaged for each group. Then \code{\link{gjack}} is called for each pair of groups. 




}
\value{
A numeric scalar of estiamted pi0.
}
%\references{
%
%
%}
\author{ Long Qu [\email{longor@iastate.edu}] }
%\note{
%
%
%}
\seealso{ 
    \code{\link{gjack}}

}
\examples{
    data(simulatedSubt)
    agjack.pi0(simulatedSubt)  ##



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric  }
%\keyword{  }
