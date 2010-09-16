\name{   pavaf1       }
\alias{    pavaf1     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
 pooling adjacent violator algorithm estimate of p-value density at 1
}
\description{
 pooling adjacent violator algorithm estimate of p-value density at 1
}
\usage{
pavaf1(p,max.bins=20,bin.method=c("max","Sturges","Scott","FD"),
                discrete=FALSE,plotit=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{p-value vector}
  \item{max.bins}{max number of bins}
  \item{bin.method}{binning method}
  \item{discrete}{logical: whether p-values are discrete.}
  \item{plotit}{logical: whether results are plotted. }
  \item{...}{Other arguments to \code{hist}}
}
\details{
This function bin the p-values and then run PAVA to estimate the minimum of its density. 




}
\value{Numeric scalar
%  \item{}{              }
%  \item{}{              }
%  \item{}{              }
  
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
\keyword{univar  }
%\keyword{  }
%\keyword{  }
