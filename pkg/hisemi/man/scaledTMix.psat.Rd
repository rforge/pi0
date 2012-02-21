\name{scaledTMix.psat}
\alias{scaledTMix.psat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fits a partially saturated model to t-statistics }
\description{
        Fits two-component mixture model to t-statistics, where each t-statistic has a different mixing proportion \code{pi0}, but all t-statistics shares a common scale factor. 
}
\usage{
scaledTMix.psat(tstat, df, upper0 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A numeric vector t-statisics }
  \item{df}{ A numeric scalar or vector of degrees of freedom }
  \item{upper0}{ A numeric scalar giving the initial upper bound to search for scale factor. It will be automatically increased if the initial bound is not appropriate. }
}
\details{
  This function fits a two-component mixture model, with a central t component with probability pi0, and a scaled central t component with scale factor greater than 1. 
  The model assumes a different pi0 for each t-statistic, but a common scale factor for all t-statistics. Maximum likelihood estimates are obtained. 
  \code{upper0} is only used as a hint of the upper bound of the scale factor. If it is too small, it will be automatically increased. 
}
\value{
  A numeric scalar, being estimated scale factor, with attributes
  \item{equiv.sd.ncp }{ A numeric scalar being the equivalent standard deviation of the noncentrality parameters, given the noncentrality parameter being nonzero. }
  \item{df }{ degrees of freedom}
  \item{fit }{ results from \code{\link{optimize}}}
  \item{n2ll}{ Negative 2 times the log likelihood}
  \item{pi0}{ A numeric vector of 0 or 1, which are the maximum likelihood estimate of \code{pi0}}
}
%\references{ }
\author{ Long Qu \email{longorognol@hotmail.com}
}
%\note{}
\seealso{ \code{\link{scaledTMix.null}}, \code{\link{scaledTMix.sat}}}
\examples{
set.seed(99927220)
(tstat=rt(5,1))
scaledTMix.psat(tstat,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
