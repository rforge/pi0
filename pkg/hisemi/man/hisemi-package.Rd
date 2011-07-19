\name{hisemi-package}
\alias{hisemi-package}
\alias{hisemi}
\docType{package}
\title{
The package includes functions for fitting hierarchical semiparametric regression model to a large number of parametric test statistics. 
}
\description{
The package includes functions for fitting hierarchical semiparametric regression model to a large number of parametric test statistics
}
\details{
\tabular{ll}{
Package: \tab hisemi\cr
Type: \tab Package\cr
Version: \tab 1.0-188\cr
Date: \tab 2011-07-18\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}

\itemize{
\item{\code{\link{penLik.EMNewton}}}{ is the major interface function to be called. }
\item{\code{\link{coef.hisemit}}} extracts the estimated parameters. 
\item{\code{\link{confint.hisemit}}} returns the confidence intervals. 
\item{\code{\link{directSum}}} computes the direct sum of matrices. 
\item{\code{\link{EMupdate}}} is the EM algorithm used in \code{\link{penLik.EMNewton}}. 
\item{\code{\link{fitted.hisemit}}} extracts the fitted values. 
\item{\code{\link{logLik.hisemit}}} returns the log likelihood. 
\item{\code{\link{NRupdate}}} is the Newton-Raphson algorithm used in \code{\link{penLik.EMNewton}}. 
\item{\code{\link{plot.hisemit}}} plots the fitted model. 
\item{\code{\link{print.hisemit}}} print summary information. 
\item{\code{\link{resid.hisemit}}} returns the residuals. 
\item{\code{\link{scaledTMix.null}}} fits the null model with common pi0. 
\item{\code{\link{scaledTMix.psat}}} fits the partially saturated model with free pi0 and common scale factor. 
\item{\code{\link{scaledTMix.sat}}} fits the completely saturated model with free pi0 and free scale factor. 
\item{\code{\link{vcov.hisemit}}} returns the sandwich variance-covariance matrix. 
}
}
\author{
Long Qu \email{longor@iastate.edu}, Dan Nettleton \email{dnett@iastate.edu}, Jack C.M. Dekkers \email{jdekkers@iastate.edu}

Maintainer: Long Qu \email{longor@iastate.edu}
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ optimize }
\seealso{
\code{\link[pi0]{pi0-package}}
}
\examples{
## simulate some fake data
G=20000
sdncp=1.3
n1=n2=5
df=n1+n2-2
set.seed(54457704)
x=runif(G,1,G)
f=function(x)sin(x*pi/1000)+1
Pi.i=1/(1+exp(f(x)))
Z.i=rbinom(G,1,1-Pi.i)
t0.i=rt(G,df)
ncp.i=rnorm(G,0,sdncp)
t1.i=rt(G,df,ncp.i)
t.i=ifelse(Z.i==0,t0.i,t1.i)

## fit model
(plfit=penLik.EMNewton(t.i, x, df, plotit=FALSE))
(plfit0=scaledTMix.null(t.i, df))

\dontrun{
plot(plfit)
plot(t.i, plfit$lfdr, pch='.')
lines(sort(t.i), plfit0$lfdr[order(t.i)], col=2, lwd=3)
}
}
