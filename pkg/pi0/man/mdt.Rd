\name{marginal.dt}
\alias{marginal.dt}
\alias{marginal.dt.parncp}
\alias{marginal.dt.nparncp}
\alias{marginal.dt.sparncp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimated arginal density of t-statistics
}
\description{
        Estimated arginal density of t-statistics from \code{ncpest} class
}
\usage{
marginal.dt(obj,...)
\method{marginal.dt}{parncp}(obj,...)
\method{marginal.dt}{nparncp}(obj, ...)
\method{marginal.dt}{sparncp}(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object of \code{ncpest} (\code{\link{nparncp}} or \code{\link{parncp}} }
%  \item{observed}{logical; if \code{TRUE}, then densities at observed \code{obj$data$tstat} are evaluated. }
  \item{\dots}{ Other argument passed to \code{\link{dtn.mix}}, most notably, the \code{approximaiton} argument }
}
\details{
When \code{obj$data$df} are all equal to each other, a single marginal density is clearly defined for all \code{obj$data$tstat}. Otherwise, 
the marginal density is defined as a discrete mixture of densities, one for each distinct degree of freedom, 
with mixing proportion based on that of \code{obj$data$df}. 
}
\value{
 A function of one argument, i.e., the marginal density function. 
}
\references{ 
}
\author{ Long Qu \email{longor@iastate.edu}
}
\note{ 
}
\seealso{ \code{\link{parncp}}. \code{\link{nparncp}} ,\code{\link{sparncp}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
