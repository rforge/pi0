\name{scaledTMix.sat}
\alias{scaledTMix.sat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fits saturated model to t-statistics }
\description{
  Fit saturated model to t-statistics, i.e., a two-component mixture model (a central t and a scaled central t with scale greater than 1) to \emph{each} t-statistics separately. 
}
\usage{
scaledTMix.sat(tstat, df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A numeric vector of t-statistics. }
  \item{df}{ A numeric scalar or vector of the same length as \code{tstat}, giving the degrees of freedom for the \code{tstat}. }
}
\details{
  This functions assumes each t-statistics coming from either a central t-distribution or a scaled central t-distribution. Each t-statistic has a different mixing proportion pi0, whose maximum likelihood estimate will be either 0 or 1. Each t-statistic has a different scale parameter. If pi0=1, the scale parameter will be 1; if pi0=0, the scale parameter will be greater than 1. 
}
\value{
  A numeric vector of estimated scale parameters, with two attributes
  \item{pi0 }{ A numeric vector of estimated pi0}
  \item{logLik }{ A numeric vector of log likelihood}
}
%\references{ }
\author{ Long Qu \email{longor@iastate.edu}
}
\note{ 
Whenever the absolute value of the \code{tstat} is less than 1, \code{pi0} will be esitmated to be 1 and the scale will also be 1. Otherwise, the \code{pi0} will be estimated to be 0 and scale will be the absolute value of \code{tstat}.
}
\seealso{ \code{\link{scaledTMix.null}}, \code{\link{scaledTMix.psat}}}
\examples{
set.seed(99927220)
(tstat=rt(5,1))
scaledTMix.sat(tstat,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
