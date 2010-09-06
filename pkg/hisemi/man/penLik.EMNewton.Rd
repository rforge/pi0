\name{penLik.EMNewton}
\alias{penLik.EMNewton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fits hierarchical semiparametric regression model to t-statistics }
\description{
  Fits hierarchical semiparametric regression model to t-statistics
}
\usage{
penLik.EMNewton(tstat, x, df, spar = 10^(-1:8), nknots = 100, starts, tuning.method = c("NIC", "CV"), 
        cv.fold = 5, optim.method = c("nlminb", "BFGS", "CG", "L-BFGS-B", "Nelder-Mead", "SANN", "NR"), 
        logistic.correction = TRUE, em.iter.max = 10, em.beta.iter.max = 1, newton.iter.max = 1500, 
        scale.conv = 0.001, lfdr.conv = 0.001, NPLL.conv = 0.001, debugging = FALSE, plot.it = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A numeric vector t-statistics }
  \item{x}{A numeric matrix of covariates, with \code{nrow(x)} being \code{length(tstat)}}
  \item{df}{ A numeric scalar or vector of degrees of freedom }
  \item{spar}{ A numeric vector of smoothing parameter lambda }
  \item{nknots}{ A numeric scalar of number of knots }
  \item{starts}{ An optional numeric vector of starting values}
  \item{tuning.method}{ Either \code{'NIC'} or \code{'CV'}, specifying the method to choose the tuning parameter \code{spar} }
  \item{cv.fold}{ A numeric scalar of the fold for cross-validation. Ignored if \code{tuning.method='NIC'}. }
  \item{optim.method}{ A character scalar specifying the method of optimization }
  \item{logistic.correction}{ A logical scalar specifying whether or not the effective number of parameters should be corrected using a logistic curve }
  \item{em.iter.max}{ A numeric scalar specifying the maximum number of EM iterations. If being \code{Inf}, then EM algorithm is used. If being \code{0}, then Newton method is used. Otherwise, EM algorithm is used initially, followed by Newton method. }
  \item{em.beta.iter.max}{ A numeric scalar specifying the maximum number of iterations in the maximization step for the beta parameters in the EM algorithm. If being \code{Inf}, the original EM is used. If being 1 or other numbers, the generalized EM algorithm is used. }
  \item{newton.iter.max}{ A numeric scalar specifying the maximum number of iterations in Newton method. }
  \item{scale.conv}{ A small numeric scalar specifying the convergence criterion for the scale parameter.  }
  \item{lfdr.conv}{ A small numeric scalar specifying the convergence criterion for the local false discovery rates. }
  \item{NPLL.conv}{ A small numeric scalar specifying the convergence criretion for the negative penalized log likelhood.  }
  \item{debugging}{ A logical scalar. If \code{TRUE}, then \code{dump.frame} will be called whenever error occurs. }
  \item{plot.it}{ A locgical scalr specifying whether a plot should be generated. }
  \item{\dots}{ Currently not used. }
}
\details{
%  ~~ If necessary, more details than the description above ~~
}
\value{
        An list of class \code{hisemit}:
\item{lfdr}{A numeric vector of local false discovery rates.}
\item{model}{A list of \code{tstat}, \code{df} and \code{x}, which are the same as arguments}
\item{scale.fact}{A list with \itemize{
        \item{scale.fact}{Scale factor}
        \item{sd.ncp}{Equivalent standard deviation of noncentrality parameters}
        \item{r}{A reparameterization of \code{scale.fact}}
        \item{t.cross}{\code{sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1))))}}}}
\item{pi0}{A numeric vector of mixing proportions for the central t component}
\item{tuning}{A list with \itemize{
        \item{mean}{Mean criterion}
        \item{var}{Variance of criterion across observations}
        \item{grp}{Cross-validation group membership}
        \item{method}{The \code{tuning.method} used.}
        \item{final}{The minimum mean criterion}}}
\item{spar}{A list with \itemize{
        \item{all}{All smoothing parameters searched}
        \item{final}{The smoothing parameter used}
        \item{final.idx}{The index of the final \code{spar}}}}
\item{enp}{A list with \itemize{
        \item{raw}{Raw effective number of parameters}
        \item{logistic}{effective number of parameters after fitting logistic curve as a correction}
        \item{final}{The effective nubmer of parameters in the final model}
        \item{good.idx}{The index of the selected effective number of parameters}}}
\item{fit}{A list with \itemize{
        \item{intercept}{The fitted intercept}
        \item{covariate.idx}{The index of covariates}
        \item{f.covariate}{Each additive smooth function evaluated at the covariates}
        \item{f}{Fitted smoothing funciton}
        \item{beta}{Estimated regression coefficients}
        \item{H}{Expanded design matrix}
        \item{asym.vcov}{Asymptotic variance-covariance matrix for estimated parameters}}}
\item{NPLL}{A list with \itemize{
        \item{NPLL}{Negative penalized log likelihood}
        \item{logLik}{Log likelihood}
        \item{penalty}{Penalty term}
        \item{saturated.ll}{Saturated log likelihood}}}
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers. A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. manuscript.  }
\author{ Long Qu \email{longor@iastate.edu} }
\note{ 
}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{resid.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{plot.hisemit}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line