\name{print.varComp}
\alias{print.varComp}
\alias{summary.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary and printing of varComp objects
}
\description{
Summary and printing functions of \code{varComp} objects. 
}
\usage{
\method{print}{varComp}(x, ...)
\method{summary}{varComp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{
	A \code{varComp} object
}
  \item{\dots}{
	Not used in \code{print.varComp}. For \code{summary.varComp}, the \code{\dots} are passed to \code{fixef.varComp}. 
}
}
\details{
\code{print.varComp} will print the fixed effect estimates, variance parameter estimates, and the number of observations. 

\code{summary.varComp} will expand the \code{fixef} component by calling \code{\link{fixef.varComp}}. 
}
\value{
\code{invisible(x)} for \code{print.varComp}. 

An object of class \code{c("summary.varComp", "varComp")} for \code{summary.varComp}. 
}
\author{
Long Qu
}
\keyword{ methods }
\keyword{ print}% __ONLY ONE__ keyword per line
