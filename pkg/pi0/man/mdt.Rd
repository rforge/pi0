\name{marginal.dt}
\alias{marginal.dt}
\alias{marginal.dt.parncpt}
\alias{marginal.dt.nparncpt}
\alias{marginal.dt.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimated arginal density of t-statistics
}
\description{
        Estimated arginal density of t-statistics from \code{ncpest} class
}
\usage{
marginal.dt(obj,...)
\method{marginal.dt}{parncpt}(obj,...)
\method{marginal.dt}{nparncpt}(obj, ...)
\method{marginal.dt}{sparncpt}(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object of \code{ncpest} (\code{\link{nparncpt}} or \code{\link{parncpt}} }
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
Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)

}
\author{ Long Qu \email{longor@iastate.edu}
}
%\note{ 
%}
\seealso{ \code{\link{parncpt}}. \code{\link{nparncpt}} ,\code{\link{sparncpt}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
