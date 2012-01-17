\name{nparncpt}
\Rdversion{1.1}
\alias{nparncpt}
\alias{nparncpt.sqp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nonparametric estimation of noncentrality parameters
}
\description{The functions use Gaussian basis functions to estimate the noncentrality parameters (ncp) from a large number of t-statistics. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nparncpt(tstat, df, ...)
nparncpt.sqp(tstat, df, penalty=c('3rd.deriv','2nd.deriv','1st.deriv'), lambdas=10^seq(-1,5,by=1), 
        starts, IC=c('BIC','CAIC','HQIC','AIC'), K=100, bounds=quantile(tstat,c(.01,.99)), 
        solver=c('solve.QP','lsei','ipop','LowRankQP'),plotit=FALSE, verbose=FALSE, approx.hess=TRUE, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ Numeric vector of noncentrality parameters
}
  \item{df}{ Numeric vector of degrees of freedom
}
  \item{penalty}{ Character vector. Penalty on either integrated squared first-order derivatives (\code{1st.deriv}) or second-order derivatives (\code{2nd.deriv}) 
  of the estimated density function of ncp. Note that only the first element is used. 
}
  \item{lambdas}{Numeric vector of smoothness tuning parameter \code{lambda} to be tried. The one that minimizes NIC will be chosen. 
}
  \item{starts}{ Optional numeric vector of starting values. If missing, \code{parncpt} will be called with \code{zeromean} set to \code{FALSE} 
  to get an initial esimate of \code{pi0}. And the starting values (\code{theta}) will be set all equal to each other and sum to \code{1-pi0}. 
  Note that this is the starting value for the largest \code{lambdas} only. For smaller \code{lambdas}, the estimates from larger \code{lambdas}
  will be used as starting values (i.e., warm start). 
}
%  \item{smooth.enp}{Additional smoothing on effective number of parameters (ENP) over \code{log10(lambdas)}. This is not well tested and should be use with care. It is rarely needed except for \emph{very} small \code{lambdas}.
%}
  \item{IC}{Character; one of \code{AIC}, \code{BIC}, \code{CAIC}, \code{HQIC}, specifying the factor multiplied to the ENP in computing Information Criterion (IC). 
}
  \item{K}{The number of basis Gaussian density functions. 
}
  \item{bounds}{A numeric vector of length 2, giving the approximate bounds where most of the probability of ncp lies. 
}
  \item{solver}{Character. The name of the function for solving quadratic programming problems. Note that \code{ipop} and \code{kernlab} are not very reliable. \code{solve.QP} is faster but \code{lsei} is more stable. 
}
  \item{plotit}{logical; indicating if \code{\link{plot.nparncpt}} should be called after estimation. This is always recommended before accepting the results. 
}
  \item{verbose}{logical; if \code{TRUE}, extensive messages will be printed. 
}
  \item{approx.hess}{either logical or a number between 0 and 1. This helps in reducing time in evaluating the hessian matrix. If it is set to \code{TRUE}, 
  for the kth Gaussian basis function and the gth \code{tstat}, the marginal t-statistic density evaluated at this \code{tstat} will be set to zero 
  if it is below the average of all \code{K*length(tstat)} such values. If it is set to FALSE or 0, then none of the density will be treated as zero, 
  no matter how small they are. If it is set to a number between 0 and 1, values below this quantile will be treated as zero. Note that this approximation 
  only affects the computation of hessian matrix, which does not need to be exact in an optimization routine. Hence, a reasonable sparseness speeds 
  up computation of a hessian matrix but might increase the number of iterations to converge. Set this to \code{TRUE} seems a reasonable trade-off between 
  the two effects and usually saves computing time. 
}
  \item{\dots}{other paramters passed to \code{\link{dtn.mix}}. Usually, the \code{approximation} argument. 
}
}
\details{ \code{nparncpt} is a wrapper for \code{nparncpt.sqp}, the latter of which uses a sequential quadratic programming algorithm to find the mixing proportions
of the basis Gaussian density functions.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A list with class attribute \code{c("nparncpt", "ncpest")}
        \item{pi0}{estimated proportion of true nulls}
        \item{mu.ncp}{mean of ncp}
        \item{sd.ncp}{SD of ncp}
        \item{logLik}{an object of class \code{logLik}. The associated \code{df} is the estimated effective number of parameters (enp). The log likelihood is also penalized likelihood. See also \code{\link{logLik.ncpest}} and \code{\link{AIC}}.}
        \item{enp}{estimated ENP}
        \item{par}{estimated parameters \code{theta}}
        \item{lambda}{the lambda that minimizes NIC}
        \item{gradiant}{analytic gradiant at the estimate}
        \item{hessian}{analytic hessian at the estimate}
        \item{beta}{estimated mixing proportions for the NCP distribution}
        \item{IC}{the information criterion specified by the user}
        \item{all.mus}{mean of each basis Gaussian density}
        \item{all.sigs}{SD of each basis Gaussian density}
        \item{data}{a list of \code{tstat} and \code{df}}
        \item{i.final}{the index of \code{lambdas} that minimizes NIC}
        \item{all.pi0s}{estimated pi0 for each lambda}
        \item{all.enps}{ENP for each lambda}
        \item{all.thetas}{parameter estimates for each lambda}
        \item{all.nics}{Network information criterion (NIC) for each lambda}
        \item{all.nic.sd}{SD of NIC for each lambda}
        \item{all.lambdas}{the \code{lambdas} argument itself}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Qu L, Nettleton D, Dekkers JCM. (2012) Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. Biometrics (in press).
%% ~put references to the literature/web site here ~
}
\author{Long Qu 
%%  ~~who you are~~
}
\note{
\code{df} could be \code{Inf} for z-tests.  When this is the case, \code{approximation} is ignored. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncpt}}, \code{\link{sparncpt}}, 
\code{\link{fitted.nparncpt}}, \code{\link{plot.nparncpt}}, \code{\link{summary.nparncpt}},
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
\keyword{ optimize}