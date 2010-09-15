\name{   ncpest.iter       }
\alias{    ncpest.iter     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
iterative call to the ncpest function
}
\description{
A wrapper to iteratively call the ncpest function

}
\usage{
ncpest.iter(p,estimates=c("all","compromise","pi0","f1"),iter=2,weights,
           eps=1e-6,keep.cdf=TRUE,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{p-value vector}
  \item{estimates}{Character: what to estimate              }
  \item{iter}{max number of iterations}
  \item{weights}{bin weights}
  \item{eps}{Small tolerance number}
  \item{keep.cdf}{Logical, as in \code{ncpest}}
  \item{\dots}{other arguments to \code{ncpest}}
}
%\details{
%
%
%
%
%
%}
\value{An \code{ncpest} object
%  \item{}{              }
%  \item{}{              }
%  \item{}{              }
  
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
\keyword{smooth  }
%\keyword{  }
%\keyword{  }
