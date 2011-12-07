\name{reflect and fold}
\alias{reflect}
\alias{fold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reflection and folding of a function about zero
}
\description{
\code{reflect} reflects the function defined on the nonnegative real line about zero to get a function defined on the whole real line, 
and then divide it by 2. \code{fold} folds a function defined on the whole real line at zero to get a function defined only on the non-negative 
real line. 
}
\usage{
reflect(f, ...)
fold(f, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ the function to be reflected or folded}
  \item{\dots}{ other arguments passed to \code{f} }
}
\details{
See examples.
}
\value{ the new function
}
%\references{ }
\author{ Long Qu }
%\note{
%}
%\seealso{  }
\examples{
## reflect function is currently defined as
    function(x,...) ifelse(x>0, f(x,...), f(-x,...))/2

## fold function is currently defined as
    function(x,...) ifelse(x>=0, f(x,...)+f(-x,...), 0)

## double exponential pdf
ddexp=reflect(dexp)

## folded normal pdf
dfnorm=fold(dnorm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
