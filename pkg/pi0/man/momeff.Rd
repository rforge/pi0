\name{   momeff       }
\alias{    momeff     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
method of moment estimate of the effect size model

}
\description{
method of moment estimate of the effect size model
}
\usage{

momeff(tstat,n1,n2=n1,zero.mean,gamma2,lower.df=6.1,upper.df=100,approx=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{t-statistics              }
  \item{n1}{Treatment 1 sample size}
  \item{n2}{Treatment 2 sample size}
  \item{zero.mean}{ logical, indicating whether a zero-mean normal model is assumed. }
  \item{gamma2}{Gamma square parameter, i.e., variance of effect sizes.}
  \item{lower.df}{  lower bound of degrees of freedom, in case of n1 is missing            }
  \item{upper.df}{  upper bound of degrees of freedom, in case of n1 is missing            }
  \item{approx}{  logical, indicating if no exact solutions are available, whether approx. solutions are returned.  }
}
%\details{
%
%
%
%
%
%}
\value{Numeric vector of parameter estimates
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
