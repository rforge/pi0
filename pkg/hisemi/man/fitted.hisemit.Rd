\name{fitted.hisemit}
\alias{fitted.hisemit}
\alias{fitted.values.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract fitted values from a hisemit object }
\description{
  Extract fitted values from a \code{hisemit} object
}
\usage{
\method{fitted}{hisemit}(object, fitted.type = c("lfdr", "fpp", "pi0", "f"), 
        gene.list, component, ...)
\method{fitted.values}{hisemit}(object, fitted.type = c("lfdr", "fpp", "pi0", "f"), 
        gene.list, component, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object }
  \item{fitted.type}{ One of \code{c("lfdr", "fpp", "pi0", "f")}. See details. }
  \item{gene.list}{ Not implemented  }
  \item{component}{ Specifying which additive component is extracted. See details.  }
  \item{\dots}{ Not used.   }
}
\details{
        For the \code{fitted.type} argument, \itemize{
        \item{\code{'lfdr'}}{ specifys local false discovery rates to be extracted.}
        \item{\code{'fpp'}}{ specifys false positive proportion to be extracted.}
        \item{\code{'pi0'}}{ specifys prior probability of null hypothesis to be extracted.}
        \item{\code{'f'}}{ specifys the smooth function to be extracted.}
}
        When \code{fitted.type='f'}, the \code{component} specifies which component of the fitted smooth function to be extracted. \cr
        If \code{component} is missing, the overal function is returned. \cr
        If \code{component='intercept'}, only the intercept term is returned. \cr
        If \code{component} is numeric, then it specifies the smooth function corresponding to which covariate to be extracted. 
}
\value{
        A numeric vector
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
\keyword{regression}
\keyword{models}
