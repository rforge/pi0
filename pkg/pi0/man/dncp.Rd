\name{dncp}
\alias{dncp}
\alias{dncp.parncpt}
\alias{dncp.nparncpt}
\alias{dncp.sparncpt}
\alias{dncp.nparncpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density of noncentrality parameters}
\description{
These functions return the density function of noncentrality parameters, from \code{ncpest} objects
}
\usage{
dncp(obj, ...)
\method{dncp}{parncpt}(obj, fold=FALSE, ...)
\method{dncp}{nparncpt}(obj, fold=FALSE, ...)
\method{dncp}{sparncpt}(obj, fold=FALSE, ...)
\method{dncp}{nparncpp}(obj, reflect=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object from \code{ncpest}, from which noncentrality parameter density to be extracted }
  \item{fold}{Logical: if \code{TRUE}, then the density of noncentrality paramters is folded about zero to give the density of absolute noncentrality parameters. 
}
  \item{reflect}{Logical: if \code{TRUE}, then the density of absolute noncentrality parameters is reflected about zero.
}
  \item{\dots}{ Further arguments.  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
        A function of one argument
}
\references{ 
Ruppert, Nettleton, Hwang. 2007. Exploring the Information in $p$-values for the Analysis and Planning of Multiple-test Experiments. Biometrics. 63. 483-495.

Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)
}
\author{ Long Qu \email{longor@iastate.edu}
%~~who you are~~ 
}
\note{ 
\code{dncp.nparncpp} is not yet implemented. 
}
\seealso{ \code{\link{parncpt}}. \code{\link{nparncpt}} ,\code{\link{sparncpt}}. \code{\link{nparncpp}}
}
%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function(obj, ...) UseMethod("dncp")
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ distribution }% __ONLY ONE__ keyword per line
