\name{logLik.hisemit}
\alias{logLik.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract the log likelihood from a hisemit object }
\description{
        Extract the log likelihood from a \code{hisemit} object.
}
\usage{
\method{logLik}{hisemit}(object, take.sum = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object }
  \item{take.sum}{ A logical scalar, indicating whether total log likelihood or the log likelihood of each data point should be extracted. }
  \item{\dots}{ Currently not used}
}
%\details{
%}
\value{
        An object of class \code{'logLik'}. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\author{ Long Qu \email{lqu@biostatsolutions.com} }
%\note{ 
%}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{resid.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{print.hisemit}},
\code{\link{penLik.EMNewton}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
