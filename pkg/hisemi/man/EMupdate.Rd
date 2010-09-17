\name{EMupdate}
\alias{EMupdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utility function performing EM algorithm updates }
\description{
Utility function performing EM algorithm updates for \code{\link{penLik.EMNewton}}
}
\usage{
EMupdate(starts, nLogLik.pen, optim.method, H, tstat, df, dt0, spar.Pen.mat, 
        em.iter.max = 10, em.beta.iter.max = 1, scale.conv = 0.001, lfdr.conv = 0.001, NPLL.conv = 0.001, debugging = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{starts}{ A numeric vector of starting values, in \code{'r'} parameterization of the \code{scale.fact}or. }
  \item{nLogLik.pen}{ A function computing negative penalized log likelihod }
  \item{optim.method}{ One of \code{c("BFGS","CG","L-BFGS-B","Nelder-Mead", "SANN", 'nlminb','NR')}, the method used for optimization. }
  \item{H}{ Design matrix }
  \item{tstat}{ A numeric vector of t-statistics}
  \item{df}{ A numeric scalar or vector of degrees of freedom }
  \item{dt0}{ A numeric vector of the central t-density evaluated at the t-statistics }
  \item{spar.Pen.mat}{ smoothing parameter times the penalty matrix }
  \item{em.iter.max}{ Maximum number of EM iterations }
  \item{em.beta.iter.max}{ Maximum number of iterations in maximization step with respect to regression coefficients }
  \item{scale.conv}{ A small numeric scalar specifying the convergence criterion for the scale parameter }
  \item{lfdr.conv}{ A small numeric scalar specifying the convergence criterion for the local false discovery rates }
  \item{NPLL.conv}{ A small numeric scalar specifying the convergence criterion for the negative penalized log likelihood }
  \item{debugging}{ A logical scalar indicating whether debugging mode of the code should be run }
}
%\details{
%}
\value{
A numeric vector of updated parameter estimates. The scale factor is in the \code{log(scale.fact-1)} parameterization. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\author{ Long Qu \email{longor@iastate.edu} }
%\note{ 
%}
\seealso{ \code{\link{penLik.EMNewton}}, \code{\link{NRupdate}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
