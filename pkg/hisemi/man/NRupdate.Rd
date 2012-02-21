\name{NRupdate}
\alias{NRupdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utility function performing Newton-Raphson algorithm updates }
\description{
Utility function performing Newton-Raphson algorithm updates for \code{\link{penLik.EMNewton}}
}
\usage{
NRupdate(f, starts, gradient, hessian, ..., ridge0 = 1e-06, 
        tolerance = sqrt(.Machine$double.eps), 
        iter.max = 1500, halving.max = Inf, relative = FALSE, 
        return.hessian = FALSE, debugging = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ Objective function to be minimized }
  \item{starts}{ A numeric vector of starting values }
  \item{gradient}{ The gradient function of \code{f} }
  \item{hessian}{ The Hessian function of \code{f} }
  \item{\dots}{ Additional arguments to be passed to \code{f} }
  \item{ridge0}{ A small ridge factor; obsolete. Current version uses \code{\link{nearPD}} to stabilize hessian }
  \item{tolerance}{ A small numeric scalar giving the convergence criterion }
  \item{iter.max}{ Maximum number of iterations }
  \item{halving.max}{ Maximum number of step-halfing }
  \item{relative}{ A logical scalar indicating if relative convergence should be checked. }
  \item{return.hessian}{ A logical scalar indicating if the final Hessian matrix is returned. }
  \item{debugging}{ A logical scalar indicating if the debuggging mode of the code should be run. }
}
%\details{
%}
\value{
    A numeric vector of updated parameters, with attributes
    \item{\code{'objective'}}{ The final evaluated objective function}
    \item{\code{'gradient'}}{ The final gradient vector}
    \item{\code{'iter'}}{ The number of iterations}
    \item{\code{'hessian'}}{ The final Hessian matrix, only available if \code{return.hessian=TRUE}. }
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers (2012) A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. Biometrics (to appear)}
\author{ Long Qu \email{longorognol@hotmail.com} }
%\note{ 
%}
\seealso{ \code{\link{penLik.EMNewton}}, \code{\link{EMupdate}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
