\name{coef.hisemit}
\alias{coef.hisemit}
\alias{coefficients.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extracts fitted parameters from a hisemit object }
\description{
        Extracts fitted parameters from a \code{hisemit} object
}
\usage{
\method{coef}{hisemit}(object, scale.parameterization = c("r", 
        "scale.factor", "sd.ncp"), ...)
\method{coefficients}{hisemit}(object, scale.parameterization = c("r", 
        "scale.factor", "sd.ncp"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object }
  \item{scale.parameterization}{ One of \code{"r"}, \code{"scale.factor"}, \code{"sd.ncp"}. See details.  }
  \item{\dots}{ Not used.  }
}
\details{
        For the scale parameter, there are three parameterizations. \itemize{
        \item{\code{'scale.factor'}}{ means the multiplicative scaling factor (greater than 1). }
        \item{\code{'sd.ncp'}}{ means the equivalent standard deviation of the noncentrality parameters. }
        \item{\code{'r'}}{ means \code{log(scale.fact-1)}, which is in the range of whole real line. }
}
}
\value{
        A numeric vector of estimated parameters
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\author{ Long Qu \email{longor@iastate.edu} }
%\note{ 
%}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{resid.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{print.hisemit}},
\code{\link{penLik.EMNewton}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{models}
