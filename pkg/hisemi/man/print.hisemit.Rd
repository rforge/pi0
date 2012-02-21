\name{print.hisemit}
\alias{print.hisemit}
\alias{summary.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print a summary of a hisemit object }
\description{
  Print summaries of a \code{hisemit} object. 
}
\usage{
\method{print}{hisemit}(x, ...)
\method{summary}{hisemit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,object}{ The \code{hisemit} object to be printed. }
  \item{\dots}{ Additional arguments to be passed to \code{\link{print.default}} or \code{\link{summary.default}} for the \code{list} object.  }
}
\details{
  Currently, the function directly calls the corresponding method for the \code{list} object
}
\value{
The same as the results from the corresponding method for the \code{list} object. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers (2012) A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. Biometrics (to appear)}
\author{ Long Qu \email{longorognol@hotmail.com} }
%\note{ 
%}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{residuals.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{print.hisemit}},
\code{\link{penLik.EMNewton}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{methods }% __ONLY ONE__ keyword per line
