\name{confint.hisemit}
\alias{confint.hisemit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Wald-type asymptotic confidence intervals from a hisemit object }
\description{
Extract Wald-type asymptotic confidence intervals from a \code{hisemit} object.
}
\usage{
\method{confint}{hisemit}(object, parm = c("lfdr", "fpp", "beta", "scale.fact", 
        "sd.ncp", "r", "coef", "pi0", "f"), level = 0.95, component, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hisemit} object. }
  \item{parm}{ One of \code{c("lfdr", "fpp", "beta", "scale.fact", "sd.ncp",} \codee{"r", "coef", "pi0", "f")}. See details. }
  \item{level}{ A numeric scalar between 0 and 1, specifying the level of confidence. }
  \item{component}{Specifying which additive component to be extracted. See details. }
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
      }
        The \code{component} specifies which component of the fitted smooth function to be extracted. \cr
        If \code{component} is missing, the overal function is returned. \cr
        If \code{component='intercept'}, only the intercept term is returned. \cr
        If \code{component} is numeric, then it specifies the smooth function corresponding to which covariate to be extracted. 
      
}
\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. 
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