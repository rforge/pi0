\name{vcov.hisemit}
\alias{vcov.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract the asymptotic variance-covariance matrix of a hisemit object }
\description{
  Extract the asymptotic variance-covariance matrix of a \code{hisemit} object
}
\usage{
\method{vcov}{hisemit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object. }
  \item{\dots}{ Currently not used.  }
}
\details{
        Variance-covariance matrix for the fitted parameters. 
}
\value{
        A numerical matrix. 
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
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line