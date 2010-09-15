\name{   cond.cdf       }
\alias{    cond.cdf     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
conditional cdf of p-values given ncp

}
\description{
conditional cdf of p-values given ncp

}
\usage{
cond.cdf(p.eval,ncp,test=c("t","z"),alternative=c("two.sided","less","greater"),
    df=if(test=="z")Inf else df,keep.cdf=TRUE,suppressWarnings=TRUE) #,save=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.eval}{numeric vector, at which the conditional CDF is evaluated.               }
  \item{ncp}{Numeric vector of noncentrality parameters              }
  \item{test}{Either \code{t}-test or \code{z}-test. }
  \item{alternative}{The same as in \code{\link{t.test}}. }
  \item{df}{The degree of freedom. }
  \item{keep.cdf}{Logical. If \code{TRUE}, the computed CDF will be stored in a hidden environment \code{.pi0cdfp} within the global environment for later re-use. }
  \item{suppressWarnings}{Logical, indicating if warnings are suppressed}
}
%\details{
%
%
%
%
%
%}
\value{

A numeric matrix  
}
\references{
Ruppert, Nettleton, Hwang. 2007. Exploring the Information in $p$-values for the Analysis and Planning of Multiple-test Experiments. Biometrics. 63. 483-495.

}
\author{ Long Qu [\email{longor@iastate.edu}] }
%\note{
%
%
%}
%\seealso{ 
%
%
%}
%\examples{
%\dontrun{
%
%
%}
%
%
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution  }
%\keyword{  }
%\keyword{  }
