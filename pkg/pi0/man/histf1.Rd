\name{   histf1       }
\alias{    histf1     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
 Histogram estimator of p-value density evaluated at 1
}
\description{
 Histogram estimator of p-value density evaluated at 1. See references. 
}
\usage{
histf1(p,max.bins=20,bin.method=c("max","nmse","bootstrap",
      "Sturges","Scott","FD"),discrete=FALSE,seq.perm=FALSE,
      nboots=200,rightBoundary=FALSE,plotit=FALSE,perm.n,perm.h,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Vector of p-values              }
  \item{max.bins}{           maximum number of bins   }
  \item{bin.method}{              binning method}
  \item{discrete}{Whether p-values are discrete}
  \item{seq.perm}{Whether p-values come from sequential permutation tests}
  \item{nboots}{bootstrap sample size}
  \item{rightBoundary}{Logical; if \code{TRUE}, then the tail mean is computed from the right boundary of the chosen bin. }
  \item{plotit}{Whether to plot the histogram}
  \item{perm.n}{\code{n} for sequential permutation tests}
  \item{perm.h}{\code{h} for sequential permutation tests}
  \item{...}{Other arguments passed to \code{hist}}
}
%\details{
%
%
%
%
%
%}
\value{A numeric scalar value of estimated p-value density at 1. 
%  \item{}{              }
%  \item{}{              }
%  \item{}{              }
  
}
\references{
Nettleton, Hwang, Caldo, Wise. 2006. Estimating the number of true null hypotheses from a histogram of $p$ values. Journal of Agricultural, Biological, and Environmental Statistics. 11. 337-356.

Bancroft and Nettleton. 2009. Estimation of False Discovery Rate Using Permutation P-values with Different Discrete Null Distributions. Iowa State University Department of Statistics Preprint Series, #2009-05.

Bancroft and Nettleton. 2009. Computationally Efficient Estimation of False Discovery Rate Using Sequential Permutation P-values. Iowa State University Department of Statistics Preprint Series, #2009-04.

Liang and Nettleton. 2012. Adaptive and dynamic adaptive procedures for false discovery rate control and estimation, Journal of the Royal Statistical Society, Series B. 74. 163-182
}
\author{ Long Qu, Kun Liang  }
%\note{
%
%
%}
\seealso{ 
\code{\link{lastbin}}, \code{\link[qvalue:qvalue]{qvalue}}

}
\examples{
set.seed(9992722)
histf1(runif(5e5)^1.5) ##  [1] 0.6762
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar  }
%\keyword{  }
%\keyword{  }
