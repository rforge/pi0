\name{penLik.EMNewton}
\alias{penLik.EMNewton}
\title{ Fits hierarchical semiparametric regression model to t-statistics }
\description{
  Fits hierarchical semiparametric regression model to t-statistics
}
\usage{
penLik.EMNewton(tstat, x, df, spar = c(10^seq(-1,8,length=30), Inf),
        nknots = n.knots(length(tstat)), starts, 
	tuning.method = c("NIC", "CV"), cv.fold = 5, pen.order=1,
	poly.degree=pen.order*2-1, optim.method =
	c("nlminb", "BFGS", "CG", "L-BFGS-B", "Nelder-Mead", "SANN", "NR"), 
        logistic.correction = TRUE, em.iter.max = 10, 
        em.beta.iter.max = 1, newton.iter.max = 1500, 
        scale.conv = 0.001, lfdr.conv = 0.001, NPLL.conv = 0.001, 
        debugging = FALSE, plotit = TRUE, ...)
}
\arguments{
  \item{tstat}{ A numeric vector t-statistics }
  \item{x}{ A numeric matrix of covariates, with \code{nrow(x)} being \code{length(tstat)}}
  \item{df}{ A numeric scalar or vector of degrees of freedom }
  \item{spar}{ A numeric vector of smoothing parameter lambda }
  \item{nknots}{ A numeric scalar of number of knots }
  \item{starts}{ An optional numeric vector of starting values}
  \item{tuning.method}{ Either \code{'NIC'} or \code{'CV'}, specifying the method to choose the tuning parameter \code{spar} }
  \item{cv.fold}{ A numeric scalar of the fold for cross-validation. Ignored if \code{tuning.method='NIC'}. }
  \item{pen.order}{ A numeric scalar of the order of derivatives of which squared integration will be used as roughness penalty.}
  \item{poly.degree}{ A numeric scalar of the degree of B-splines. }
  \item{optim.method}{ A character scalar specifying the method of optimization. }
  \item{logistic.correction}{ A logical scalar specifying whether or not the effective number of parameters should be corrected using a logistic curve }
  \item{em.iter.max}{ A numeric scalar specifying the maximum number of EM iterations. If being \code{Inf}, then EM algorithm is used. If being \code{0}, then Newton method is used. Otherwise, EM algorithm is used initially, followed by Newton method. }
  \item{em.beta.iter.max}{ A numeric scalar specifying the maximum number of iterations in the maximization step for the beta parameters in the EM algorithm. If being \code{Inf}, the original EM is used. If being 1 or other numbers, the generalized EM algorithm is used. }
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
        \item{\code{scale.fact}:}{ Scale factor}
        \item{\code{sd.ncp}:}{ Equivalent standard deviation of noncentrality parameters}
        \item{\code{r}:}{ A reparameterization of \code{scale.fact}}
        \item{\code{t.cross}:}{ \code{sqrt(df*(s^(2/(df+1))-1)/(1-s^(-2*df/(df+1))))}}} where \code{s} is the \code{scale.fact}}
\item{pi0:}{ A numeric vector of mixing proportions for the central t component}
\item{tuning:}{ A list with \itemize{
        \item{\code{mean}:}{ Mean criterion}
        \item{\code{var}:}{ Variance of criterion across observations}
        \item{\code{grp}:}{ Cross-validation group membership}
        \item{\code{method}:}{ The \code{tuning.method} used.}
        \item{\code{final}:}{ The minimum mean criterion}}}
\item{spar:}{ A list with \itemize{
        \item{\code{all}:}{ All smoothing parameters searched}
        \item{\code{final}:}{ The smoothing parameter used}
        \item{\code{final.idx}:}{ The index of the final \code{spar}}}}
\item{enp:}{ A list with \itemize{
        \item{\code{raw}:}{ Raw effective number of parameters}
        \item{\code{logistic}:}{ Effective number of parameters after fitting logistic curve as a correction}
        \item{\code{final}:}{ The effective nubmer of parameters in the final model}
        \item{\code{good.idx}:}{ The index of the selected effective number of parameters}}}
\item{fit:}{ A list with \itemize{
        \item{\code{intercept}:}{ The fitted intercept}
        \item{\code{covariate.idx}:}{ The index of covariates}
        \item{\code{f.covariate}:}{ Each additive smooth function evaluated at the covariates}
        \item{\code{f}:}{ Fitted smoothing funciton}
        \item{\code{beta}:}{ Estimated regression coefficients}
        \item{\code{H}:}{ Expanded design matrix}
        \item{\code{asym.vcov}:}{ Asymptotic variance-covariance matrix for estimated parameters}}}
\item{NPLL:}{ A list with \itemize{
        \item{\code{NPLL}:}{ Negative penalized log likelihood}
        \item{\code{logLik}:}{ Log likelihood}
        \item{\code{penalty}:}{ Penalty term}
        \item{\code{saturated.ll}:}{ Saturated log likelihood}}}
}
\note{
When spar is too small, the results need to be treated cautiously. It is advisable to plot the results as a check. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers (2012) A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. Biometrics, 68(4):1168-1177}
\author{ Long Qu \email{long.qu@wright.edu} }
%\note{ 
%}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, 
\code{\link{vcov.hisemit}}, \code{\link{residuals.hisemit}}, \code{\link{logLik.hisemit}}, 
\code{\link{confint.hisemit}}, \code{\link{plot.hisemit}}, 
\code{\link{hisemi-package}}, \code{\link[pi0]{pi0-package}}
}
\examples{
%\dontrun{
# See the examples for the hisemi-package.
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
