\name{sparncpt}
\alias{sparncpt}
\alias{sparncpt.parncpt}
\alias{sparncpt.nparncpt}
\alias{sparncpt.numeric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Semiparametric density estimation for noncentrality parameters
}
\description{
Semiparametric density estimation for noncentrality parameters using the combination method of Olkin and Spiegelman (1987), 
based on fits from both \code{\link{parncpt}} and \code{\link{nparncpt}}.
}
\usage{
sparncpt(obj1, obj2, ...)
\method{sparncpt}{parncpt}(obj1, obj2, ...)
\method{sparncpt}{nparncpt}(obj1, obj2, ...)
\method{sparncpt}{numeric}(obj1, obj2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1, obj2}{ Case 1: \code{obj1} and \code{obj2} are of class \code{\link{parncpt}} and \code{\link{nparncpt}} respectively; 
  or vice versa; Case 2: \code{obj1} is a numeric vector of t-statistics and \code{obj2} is a vector degrees of freedom}
  \item{\dots}{ other arguments passed to \code{\link{dtn.mix}}, most notably the \code{approximation} argument. }
}
\details{
This is a two-component mixture of a parametric fit from \code{\link{parncpt}} and a nonparametric fit from \code{\link{nparncpt}}, 
with mixing proportion rho. If \code{obj1} and \code{obj2} are t-statistics and degrees of freedom respectively, calls to each of \code{\link{parncpt}} and \code{\link{nparncpt}}
are made and their results are used in combination.
}
\value{ a list with class \code{c('sparncpt','ncpest')}:
        \item{pi0}{estimated proportion of true nulls}
        \item{mu.ncp}{mean of ncp}
        \item{sd.ncp}{SD of ncp}
        \item{logLik}{an object of class \code{logLik}. The associated \code{df} is the estimated effective number of parameters (enp). The log likelihood is also penalized likelihood. See also \code{\link{logLik.ncpest}} and \code{\link{AIC}}.}
        \item{enp}{estimated ENP}
        \item{par}{estimated mixing proportion \code{rho}}
        \item{gradiant}{analytic gradiant at the estimate (not implemented) }
        \item{hessian}{analytic hessian at the estimate (not implemented) }
        \item{parfit}{ the fitted \code{\link{parncpt}} object}
        \item{nparfit}{the fitted \code{\link{nparncpt}} object}
}
\references{ 
I. Olkin and C. H. Spiegelman.  1987. A Semiparametric Approach to Density Estimation. Journal of the American Statistical Association. 82,399,858-865

Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)
}
\author{ Long Qu \email{lqu@biostatsolutions.com}}
%\note{ 
%}
\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, 
\code{\link{fitted.sparncpt}}, \code{\link{plot.sparncpt}}, \code{\link{summary.sparncpt}},
\code{\link{coef.ncpest}}, \code{\link{logLik.ncpest}}, \code{\link{vcov.ncpest}},
\code{\link{AIC}}, \code{\link{dncp}}
}
\examples{
\dontrun{
data(simulatedTstat)
(npfit=nparncpt(tstat=simulatedTstat, df=8)); 
(pfit=parncpt(tstat=simulatedTstat, df=8, zeromean=FALSE)); plot(pfit)
(pfit0=parncpt(tstat=simulatedTstat, df=8, zeromean=TRUE)); plot(pfit0)
(spfit=sparncpt(npfit,pfit)); plot(spfit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ smooth }% __ONLY ONE__ keyword per line
