\name{dncp}
\alias{dncp}
\alias{dncp.parncpt}
\alias{dncp.nparncpt}
\alias{dncp.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density of noncentrality parameters}
\description{
These functions return the density function of noncentrality parameters, from \code{ncpest} objects
}
\usage{
dncp(obj, ...)
\method{dncp}{parncpt}(obj, ...)
\method{dncp}{nparncpt}(obj, ...)
\method{dncp}{sparncpt}(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object from \code{ncpest}, from which noncentrality parameter density to be extracted }
  \item{\dots}{ currently not used }
}
\details{
%  ~~ If necessary, more details than the description above ~~
}
\value{
        A function of one argument
}
\references{ 
%~put references to the literature/web site here ~ 
}
\author{ Long Qu \email{longor@iastate.edu}
%~~who you are~~ 
}
\note{ 
%~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ \code{\link{parncpt}}. \code{\link{nparncpt}} ,\code{\link{sparncpt}}
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(obj, ...) UseMethod("dncp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
