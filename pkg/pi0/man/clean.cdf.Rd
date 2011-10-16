\name{   clean.cdf, save.cdf and load.cdf       }
\alias{    clean.cdf     }
\alias{    save.cdf     }
\alias{    load.cdf     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 

Saving, loading or removing pre-computed conditional cdf objects from memory or disk
}
\description{
Saving, loading or removing pre-computed conditional cdf objects from memory or disk

}
\usage{
save.cdf(path=getwd())
load.cdf(path=getwd())
clean.cdf(where=c("memory","disk","both"),path=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{where}{Clean either \code{memory} or \code{disk} or \code{both} cdf objects}
  \item{path}{Character: path where pre-computed conditional cdf RData files are stored}
}
\details{
The filename on the disk is \code{.pi0cdfp.RData}; the environment name is \code{.pi0cdfp} which contains the pre-computed conditional cdf objects. The environment itself is stored in global environment. 




}
\value{
invisible logical status for \code{save.cdf} and \code{clean.cdf}. For \code{load.cdf}, a character string from \code{load} or a \code{try-error} object.
}
\references{
Ruppert, Nettleton, Hwang. 2007. Exploring the Information in $p$-values for the Analysis and Planning of Multiple-test Experiments. Biometrics. 63. 483-495.
}
\author{ Long Qu [\email{lqu@biostatsolutions.com}] }
%\note{
%
%
%}
\seealso{ 

\code{\link{cond.cdf}}
}
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
\keyword{utilities  }
%\keyword{  }
%\keyword{  }
