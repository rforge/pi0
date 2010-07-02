\name{resid.hisemit}
\alias{resid.hisemit}
\alias{residuals.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract residuals from a hisemit object }
\description{
  Extract residuals from a \code{hisemit} object
}
\usage{
\method{resid}{hisemit}(object, residual.type = "deviance", ...)
\method{residuals}{hisemit}(object, residual.type = "deviance", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object }
  \item{residual.type}{ A character scalar specifying the type of residuals to be extracted. Currently only \code{'deviance'} is supported. }
  \item{\dots}{ Currently not used. }
}
\details{
  
}
\value{
        A numeric vector residuals. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\author{ Long Qu \email{longor@iastate.edu} }
\note{ 
}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{resid.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{print.hisemit}},
\code{\link{penLik.EMNewton}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{models}
