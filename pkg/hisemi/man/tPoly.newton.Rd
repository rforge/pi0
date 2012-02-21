\name{tPoly.newton}
\alias{tPoly.newton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fits hierarchical global polynomial regression model to t-statistics }
\description{
  Fits hierarchical global polynomial regression model to t-statistics through Newtonian algorithms. 
}
\usage{
tPoly.newton(tstat, x, df, starts, 
        pen.order=1,
        optim.method = c("nlminb", "BFGS", "CG", 
        "L-BFGS-B", "Nelder-Mead", "SANN", "NR"), 
        newton.iter.max = 1500, 
        scale.conv = 0.001, lfdr.conv = 0.001, NPLL.conv = 0.001, 
        debugging = FALSE, plotit = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A numeric vector t-statistics }
  \item{x}{ A numeric matrix of covariates, with \code{nrow(x)} being \code{length(tstat)}}
  \item{df}{ A numeric scalar or vector of degrees of freedom }
%  \item{spar}{ A numeric vector of smoothing parameter lambda }
%  \item{nknots}{ A numeric scalar of number of knots }
  \item{starts}{ An optional numeric vector of starting values. The first element is the \code{r}, i.e. \code{log(scale-1)}. The second parameter is the intercept. The remaining elements are the starting values for the B-spline coefficients (removing the first basis) for each x. When thist argument is not provided, the code starts with a global constant model that is easiest to fit, and then increase the order gradually using the warm starts from lower order fits.}
%  \item{tuning.method}{ Either \code{'NIC'} or \code{'CV'}, specifying the method to choose the tuning parameter \code{spar} }
%  \item{cv.fold}{ A numeric scalar of the fold for cross-validation. Ignored if \code{tuning.method='NIC'}. }
  \item{pen.order}{ A numeric scalar of the order of derivatives of which squared integration will be used as roughness penalty. Note: The final order of the global polynomial is always \code{pen.order-1}.}
%  \item{poly.degree}{ A numeric scalar of the degree of B-splines. }
  \item{optim.method}{ A character scalar specifying the method of optimization. }
%  \item{logistic.correction}{ A logical scalar specifying whether or not the effective number of parameters should be corrected using a logistic curve }
%  \item{em.iter.max}{ A numeric scalar specifying the maximum number of EM iterations. If being \code{Inf}, then EM algorithm is used. If being \code{0}, then Newton method is used. Otherwise, EM algorithm is used initially, followed by Newton method. }
%  \item{em.beta.iter.max}{ A numeric scalar specifying the maximum number of iterations in the maximization step for the beta parameters in the EM algorithm. If being \code{Inf}, the original EM is used. If being 1 or other numbers, the generalized EM algorithm is used. }
  \item{newton.iter.max}{ A numeric scalar specifying the maximum number of iterations in Newton method. }
  \item{scale.conv}{ A small numeric scalar specifying the convergence criterion for the scale parameter.  }
  \item{lfdr.conv}{ A small numeric scalar specifying the convergence criterion for the local false discovery rates. }
  \item{NPLL.conv}{ A small numeric scalar specifying the convergence criretion for the negative penalized log likelhood.  }
  \item{debugging}{ A logical scalar. If \code{TRUE}, then \code{dump.frame} will be called whenever error occurs. }
  \item{plotit}{ A locgical scalr specifying whether a plot should be generated. }
  \item{\dots}{ Currently not used. }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
        An list of class \code{hisemit}:
\item{lfdr:}{ A numeric vector of local false discovery rates.}
\item{model}{A  list of \code{tstat}, \code{df} and \code{x}, which are the same as arguments}
\item{scale.fact:}{ A list with \itemize{
        \item{scale.fact:}{ Scale factor}
        \item{sd.ncp:}{ Equivalent standard deviation of noncentrality parameters}
        \item{r:}{ A reparameterization of \code{scale.fact}}
        \item{t.cross:}{ \code{sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1))))}}}}
\item{pi0:}{ A numeric vector of mixing proportions for the central t component}
\item{tuning:}{ A list with \itemize{
        \item{mean:}{ Mean criterion}
        \item{var:}{ Variance of criterion across observations}
        \item{grp:}{ Cross-validation group membership}
        \item{method:}{ The \code{tuning.method} used.}
        \item{final:}{ The minimum mean criterion}}}
\item{spar:}{ A list with \itemize{
        \item{all:}{ All smoothing parameters searched}
        \item{final:}{ The smoothing parameter used}
        \item{final.idx:}{ The index of the final \code{spar}}}}
\item{enp:}{ A list with \itemize{
        \item{raw:}{ Raw effective number of parameters}
        \item{logistic:}{ Effective number of parameters after fitting logistic curve as a correction}
        \item{final:}{ The effective nubmer of parameters in the final model}
        \item{good.idx:}{ The index of the selected effective number of parameters}}}
\item{fit:}{ A list with \itemize{
        \item{intercept:}{ The fitted intercept}
        \item{covariate.idx:}{ The index of covariates}
        \item{f.covariate:}{ Each additive smooth function evaluated at the covariates}
        \item{f:}{ Fitted smoothing funciton}
        \item{beta:}{ Estimated regression coefficients}
        \item{H:}{ Expanded design matrix}
        \item{asym.vcov:}{ Asymptotic variance-covariance matrix for estimated parameters}}}
\item{NPLL:}{ A list with \itemize{
        \item{NPLL:}{ Negative penalized log likelihood}
        \item{logLik:}{ Log likelihood}
        \item{penalty:}{ Penalty term}
        \item{saturated.ll:}{ Saturated log likelihood}}}
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers (2012) A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. Biometrics (to appear)}
\author{ Long Qu \email{longorognol@hotmail.com} }
%\note{ 
%}
\seealso{ \code{\link{penLike.EMNewton}}, \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, 
\code{\link{vcov.hisemit}}, \code{\link{residuals.hisemit}}, \code{\link{logLik.hisemit}}, 
\code{\link{confint.hisemit}}, \code{\link{plot.hisemit}} 
\code{\link{hisemi-package}}, \code{\link[pi0]{pi0-package}}
}
\examples{
\dontrun{
# See the example for the package.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
