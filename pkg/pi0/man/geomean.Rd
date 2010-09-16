\name{   geomean and harmean      }
\alias{    geomean     }
\alias{    harmean     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Geometric mean and harmonic mean functions

}
\description{
Geometric mean and harmonic mean functions

}
\usage{
geomean(x)
harmean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector}
%  \item{}{              }
%  \item{}{              }
}
%\details{
%
%
%
%
%
%}
\value{numeric scalar. For \code{geomean}, this is the geometric mean of \code{x}; for \code{harmean}, this is the harmonic mean of \code{x}. 
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
\examples{
geomean(10^(1:10)) #[1] 316227.8
harmean(10^(1:10)) #[1] 90
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar  }
%\keyword{  }
%\keyword{  }
