\name{confint.hisemit}
\alias{confint.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Wald-type asymptotic confidence intervals from a hisemit object }
\description{
Extract Wald-type asymptotic confidence intervals from a \code{hisemit} object.
}
\usage{
\method{confint}{hisemit}(object, parm = c("lfdr", "fpp", "beta", "scale.fact", "sd.ncp", "r", "coef", "pi0", "f"), level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object. }
  \item{parm}{ One of \code{c("lfdr", "fpp", "beta", "scale.fact", "sd.ncp", "r", "coef", "pi0", "f")}. See details. }
  \item{level}{ A numeric scalar between 0 and 1, specifying the level of confidence. }
  \item{\dots}{ Currently not used. }
}
\details{
        The \code{parm} could be: \itemize{
        \item{\code{'lfdr'}}{ Local false discovery rates}
        \item{\code{'fpp'}}{ False positive proportions; not implemented yet. }
        \item{\code{'beta'}}{ Regression coefficients}
        \item{\code{'scale.fact'}}{ Multiplicative scale factor}
        \item{\code{'r'}}{ \code{log(scale.fact-1)} }
        \item{\code{'sd.ncp'}}{ Equivalent standard deviation of noncentrality parameters }
        \item{\code{'coef'}}{ All parameters }
        \item{\code{'pi0'}}{ Prior probability of true null hypotheses }
        \item{\code{'f'}}{ Underlying smooth function}
}}
\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. 
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
\keyword{ models }