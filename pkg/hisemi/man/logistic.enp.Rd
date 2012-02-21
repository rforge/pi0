\name{logistic.enp}
\alias{logistic.enp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a logistic curve to the raw effective number of parameters over log smoothing parameter }
\description{
        Fit a logistic curve to the raw effective number of parameters over log smoothing parameter
}
\usage{
logistic.enp(log.spar, enps, maximum, minimum = 2, eps = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{log.spar}{ A numeric vector of log smoothing parameters }
  \item{enps}{ A numeric vector of raw effective number of parameters }
  \item{maximum}{ A numeric scalar of upper bound of the effective number of parameters }
  \item{minimum}{ A numeric scalar of lower bound of the effective enumber of parameters }
  \item{eps}{ A small numeric scalar of error of tolerence }
}
\details{
The purpose is of this function is to reduce some estimation variability of effective number of parameters. The function heuristically fit a smooth logistic curve to the raw effective number of parameters. 
}
\value{
A numeric vector of fitted effective number of parameters, with attributes
\item{log.spar}{ log of smoothing parameter}
\item{rate}{ rate parameter of the logistic curve}
\item{mdis}{ the center parameter of the logistic curve}
\item{pow}{ the exponent parameter of the logistic curve}
\item{fit}{ a nonlinear least squared fit object}
\item{goodenp.idx}{ Indices of enps whose values are reasonably estimated}
\item{mode}{ The model of enps from a isotonic fit}
}
%\references{ }
\author{ Long Qu \email{longorognol@hotmail.com} }
%\note{ 
%}
\seealso{ \code{\link{nls}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
