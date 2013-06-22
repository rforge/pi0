\name{varComp}
\alias{varComp}
\alias{varComp.fit}
\alias{doFit.varComp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting variance component models
}
\description{
\code{varComp} and \code{varComp.fit} fit linear mixed-effect models where the marginal variance-covariance matrix is linear in known positive semidefinite matrices. \code{varComp} uses the usual formula interface, whereas \code{varComp.fit} is the underlying working horse. \code{doFit.varComp} performs model fitting if the object has been previously created by setting \code{doFit = FALSE} when calling \code{varComp}. 
}
\usage{
varComp(fixed, data, random, varcov, weights, subset, 
family = stats::gaussian('identity'), na.action, offset, 
control = varComp.control(...), doFit = TRUE, 
normalizeTrace = TRUE, contrasts = NULL, 
model = TRUE, X = TRUE, Y = TRUE, K = TRUE, ...)

varComp.fit(Y, X = matrix(0, length(Y), 0L), K, control = varComp.control())

doFit.varComp(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{A two-sided formula specifying the fixed effects of the model. This is the same as in \code{\link[stats:lm]{stats::lm}} or \code{\link[nlme:lme]{nlme::lme}}. 
}
  \item{data}{an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are taken from \code{environment(fixed)}, typically the environment from which \code{varComp} is called. This is the same as the \code{data} argument of \code{\link[stats:lm]{stats::lm}} or \code{\link[nlme:lme]{nlme::lme}}.
}
  \item{random}{An optional one-sided formula specifying additive random effects, of the form \code{~ z1 + z2 + z3}. Interactions are allowed. The error variance component should not be included. A warning will be given when \code{environment(fixed)} and \code{environment(random)} do not match. See details.
}
  \item{varcov}{An optional list of symmetric positive semidefinite matrices. If \code{random} is non-missing, these matrices represent the correlation matrix of each random effect. Thus the number of the random effects must match the length of \code{varcov}. If \code{random} is missing, the weighted sum of these matrices represent the contribution of random effects to the marginal variance of the response variable, with unknown weights representing variance components. 
}
  \item{weights}{
	An optional nonnegative vector of the same length as the response variable specified in \code{fixed}. When it is given, it is inversely proportional to the error variances not captured by \code{random} and \code{varcov}. This is similar to the \code{weights} argument in \code{\link[stats:lm]{stats::lm}}. 
}
  \item{subset}{An optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{family}{The same as the \code{family} argument of \code{\link[stats:glm]{stats::glm}}. However, only \code{gaussian('identity')} is supported currently. 
}
  \item{na.action}{The same as in \code{\link[stats:lm]{stats::lm}}.
}
  \item{offset}{The same as in \code{\link[stats:glm]{stats::glm}}. These offsets are assumed as known fixed effects. 
}
  \item{control}{An object from calling \code{\link{varComp.control}}. 
}
  \item{doFit}{A logical scalar, indicating whether model fitting should be performed. 
}
  \item{normalizeTrace}{A logical scalar, indicating whether the individual variance-covariance matrices should be normalized such that variance components are on the same scale. 
}
  \item{contrasts}{The same as in \code{\link[stats:lm]{stats::lm}}. 
}
  \item{model}{A logical scalar, indicating whether the model frame will be included in the result. 
}
  \item{X}{For \code{varComp}, this is a logical scalar, indicating whether the fixed-effect design matrix should be included in the result. For \code{varComp.fit} this is the optional numeric fixed effect design matrix for the model. If \code{X} is missing or a matrix with zero columns, it is assumed that \code{Y} has zero mean. 
}
  \item{Y}{For \code{varComp}, this is a logical scalar, indicating whether the response variable should be included in the result. For \code{varComp.fit}, this is a numeric vector of response variables. 
}
  \item{K}{For \code{varComp}, this is a logical scalar, indicating whether the list of variance-covariance matrices should be included in the result. These matrices are computed by pre-and-post multiplying the design matrices of random effects specified by \code{random} and the corresponding matrices specified by \code{varcov}.  For \code{varComp.fit}, this is a list of variance-covariance matrices. 
}
  \item{\dots}{When \code{control} is given, this is ignored. Otherwise, these arguments are passed to \code{\link{varComp.control}} and the result will be used as the \code{control} argument. 
}
	\item{object}{An object of class \code{varComp}. 
}

}
\details{
The variance component model is of form 
\deqn{\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{e}}{Y=X \beta + e}
where \eqn{\mathbf{e}}{e} is multivariate normally distributed with mean zero and variance-covariance matrix \eqn{\mathbf{V}}{V} being
\deqn{\mathbf{V} = \sum_{j=1}^R \sigma_j^2 \mathbf{K}_j + \sigma_e^2 \mathbf{W}}{V= \sigma{j=1}^R \sigma_j^2 K_j + \sigma_e^2 W}
in which \eqn{\mathbf{K}_j}{K_j} are known positive semidefinite matrices and \eqn{\mathbf{W}}{W} is a known diagonal positive definite matrix. In the case of random-effect modeling, the \eqn{\mathbf{K}}{K} matrices are further given by 
\deqn{\mathbf{K}_j = \mathbf{Z}_j \mathbf{G}_j \mathbf{Z}_j^\top}{K_j = Z_j G_j Z_j^T}
where \eqn{\mathbf{Z}_j}{Z_j} is the design matrix for the \eqn{j}th random-effect factor and \eqn{\mathbf{G}_j}{G_j} the known variance-covariance matrix for the \eqn{j}th random-effect vector. 

In the \code{varComp} formula interface, the \eqn{\mathbf{X}}{X} matrix and response variable are specified by the \code{fixed} argument. The optional \code{random} argument specifies the \eqn{\mathbf{Z}}{Z} matrices. When \code{random} is missing, they are assumed to be identity matrices. The \code{varcov} argument specifies the \eqn{\mathbf{G}}{G} matrices. The \code{weights} argument specifies the \eqn{\mathbf{W}}{W} matrix. 

Note that in \code{random}, kernel functions are allowed. For example, \code{random = ~ ibs(SNP)} is allowed to specify the similarities contributed by a matrix \code{SNP} through the identity-by-descent (\code{\link{IBS}}) kernel. 

Fixed-by-random interactions are allowed in \code{random}. If \code{F} is a fixed factor with 2 levels, and \code{F} has already appeared in the \code{fixed} argument, then \code{random = ~ ibs(SNP) + F:ibs(SNP)} will specify two random effects. The first contributes to the marginal variance-covariance matrix through \code{IBS(SNP)}; the other contributes to the marginal variance-covariance matrix through \code{tcrossprod(X_F2) * IBS(SNP)},  where \code{X_F2} is the 2nd column of the fixed-effect design matrix of factor \code{F} under the sum-to-zero contrast. Note that sum-to-zero contrast will be used even if another contrast is used for \code{F} in the \code{fixed} argument. This behavior corresponds to what we usually mean by fixed-by-random interaction. If one indeed needs to avoid sum-to-zero contrast in this case, it is only possible to do so by providing \code{varcov} directly without specifying \code{random}.  

Additionally, the fixed-by-random interactions in \code{random} specified by the ``\code{:}'' symbol are actually treated as if they are specified by ``\code{*}''. In other words, \code{random = ~ ibs(SNP) + F:ibs(SNP)}, \code{random = ~ F:ibs(SNP)}, and \code{random = ~ F*ibs(SNP)} are actually treated as equivalent when \code{F} has already appeared on the right-side of \code{fixed}. Again this corresponds to the usual notion of fixed-by-random interaction. If such behavior is not wanted, one can directly give the intended \code{varcov} and leave \code{random} as missing. 

Finally, intercept is not included in the \code{random} interface.

The model fitting process uses profiled restricted maximum likelihood (PREML), where the error variance is always profiled out of the REML likelihood. The non-negativity constraints of variance components are always imposed. See \code{\link{varComp.control}} for arguments controlling the modeling fitting. 
}
\value{The value of any of the three functions is a list with class \code{varComp}. Depending on whether \code{doFit} is \code{TRUE} or \code{FALSE}, the components in the list might differ. 
In either case, the result from \code{varComp} always contains the following elements: 
\itemize{
	\item \code{na.action}: the \code{na.action} used in the model frame. 
	\item \code{offset}: the same as input. 
	\item \code{contrasts}: the contrast used in the fixed-effect design matrix. 
	\item \code{xzlevels}: the levels of both fixed-effect and random-effect factors. 
	\item \code{terms}: both fixed-effect and random-effect terms. 

	\item \code{call}: the actual call. 
	\item \code{nobs}: the number of observations. 
	\item \code{control}: the same as input. 
	\item \code{random.labels}: the labels used to differentiate random effects. Note that the error variance is not included here. It is safe to check the number of variance components specified by the model by checking the length of \code{random.labels}. 
	\item \code{doFit}: The value of \code{doFit} is either \code{TRUE} or a call. If the \code{doFit} argument of \code{varComp} is set to \code{FALSE}, this component will be a call that will do the model fitting when being evaluated. If the model fitting has already been done, this component will always be \code{TRUE}.
}

In the case where \code{doFit} is \code{TRUE}, the following components will also be present in the result of \code{varComp} and \code{doFit}: 
\itemize{
	\item {\code{parms}}: a named numeric vector of parameter estimates. These are the estimated ratio of each variance component relative to the error variance. 
	\item {\code{gradients}}: a named numeric vector of gradient of the PREML function at the final parameter estimate. Because of non-negativity constraints, the gradients are not necessarily all zero at convergence. 
	\item {\code{hessian}}: a numeric matrix of Hessian matrix at parameter estimate. This is always computed through the observed information, no matter how the \code{information} argument of \code{\link{varComp.control}} is set. 
	\item {\code{sigma2}}: the estimated error variance. 
	\item {\code{varComps}}: a named numeric vector of variance components. These are the same as \code{parms} times \code{sigma2}. 
	\item {\code{n.iter}}: the number of iterations during fitting. 
	\item {\code{PREML}}: the final maximized PREML function value. 
	\item {\code{X.Q2}}: a matrix, when left-multiplied to the response variable, produces the residual contrasts. 
	\item {\code{residual.contrast}}: a numeric vector of residual contrasts with \code{length{Y} - qr(X)$rank} elements. 
	\item {\code{working.cor}}: the list of individual variance-covariance matrices, whose weighted sum is the marginal variance-covariance matrix of the residual contrast. 
}

If \code{varComp.fit} is called directly, then \code{doFit} is always \code{TRUE}, but the result will not contain any of \code{na.action}, \code{offset}, \code{contrasts}, \code{xzlevels}, and \code{terms}, as they do not apply. 

For the \code{varComp} interface, if any of \code{model}, \code{X}, \code{Y} or \code{K} is \code{TRUE}, the corresponding part will be included in the result. If \code{weights} is non-missing, it will be included in the result. 
}

\note{
	Currently, available S3 methods for the \code{varComp} class are: 
\itemize{
	\item {\code{coef}}: retrieving fixed-effect estimates, variance component estimates, or estimates of the ratios of variance components to the error variance. 
	\item {\code{model.matrix}}: retrieving the fixed-effect design matrix, equivalent random-effect design matrix, or the list of individual variance-covariance matrices. 
	\item {\code{vcov}}: computing variance-covariance matrix of fixed-effect parameters, variance parameters, or the response variable. 
	\item {\code{fixef}}: retrieving fixed-effect parameter estimates with standard errors and default hypothesis tests against zero. See also \code{\link{satterth.varComp}} and \code{\link{KR.varComp}} for computing approximate denominator of degrees of freedom for the \eqn{F}-statistics.
	\item {\code{anova}}: similar to \code{fixef}, but returning a data frame. 
	\item {\code{varComp.test}}: testing for nullity of random components (other than the error variance). 
	\item {\code{logLik}}: retrieving the maximized PREML value with degrees of freedom. 
	\item {\code{print}}: showing the object. 
	\item {\code{summary}}: summarize the object by providing standard errors and default tests against zero. 
}
}


\references{
	Qu L, Guennel T, Marshall SL. (2013+) Linear Score Tests for Variance Components in Linear Mixed Models and Applications to Genetic Association Studies. \emph{Biometrics}, under revision. 
}
\author{
Long Qu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{varComp.test}}, \code{\link{coef.varComp}}, \code{\link{model.matrix.varComp}}, \code{\link{vcov.varComp}}, \code{\link{fixef.varComp}}, \code{\link{satterth.varComp}}, \code{\link{KR.varComp}}, \code{\link{logLik.varComp}}, \code{\link{print.varComp}}, \code{\link{summary.varComp}}, \code{\link[nlme:lme]{nlme::lme}}, \code{\link[stats:lm]{stats::lm}}
}
\examples{
### Oxide/Semiconductor data example
library(nlme)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
VarCorr(lmef)
coef(vcf, 'varComp') ## same values as above
k0=tcrossprod(model.matrix(~0+Lot,Oxide))
k1=tcrossprod(Oxide$Source==1)*k0
k2=tcrossprod(Oxide$Source==2)*k0
k3=tcrossprod(model.matrix(~0+Lot:Wafer, Oxide))
## unequal variance across Source for Lot effects, in a preferred parameterization:
(vcf1 = varComp(Thickness~Source, Oxide, varcov=list(S1Lot=k1, S2Lot=k2, `Lot:Wafer`=k3))) 
## unequal variance across Source for Lot effects, in a different parameterization:
(vcf2 = varComp(Thickness~Source, Oxide, varcov=list(Lot=k0, S2Lot=k2, `Lot:Wafer`=k3)))
## unequal variance across Source for Lot effects, but in a poor parameterization that 
##   turns out to be the same as vcf after fitting.
(vcf3 = varComp(Thickness~Source, Oxide, varcov=list(Lot=k0, S1Lot=k1, `Lot:Wafer`=k3)))  
logLik(vcf)
logLik(vcf1)
logLik(vcf2)  ## the same as vcf1
logLik(vcf3)  ## the same as vcf
## fixef-effect only
vcf0 = varComp(Thickness~Source, Oxide)
summary(vcf0)
summary(lmf<-lm(Thickness~Source, Oxide))
vcf00 = varComp(Thickness~0, Oxide)
summary(vcf00)
summary(lmf0<-lm(Thickness~0, Oxide))



### Genetics example	
trt=gl(2, 15)
set.seed(2340)
dat=data.frame(trt=trt)
dat$SNP=matrix(sample(0:2, 120, replace=TRUE), 30)
dat$Y  = as.numeric(trt)+rnorm(30) + dat$SNP\%*\%rnorm(4)
(vcf0 = varComp(Y~trt, dat, ~ibs(SNP)))
(vcf00 = varComp(Y~trt, dat, varcov = list(`ibs(SNP)`=IBS(dat$SNP)))) ## same as above
(vcf1 = varComp(Y~trt, dat, ~ibs(SNP):trt)) ## two variance components

dat$trt[1]=NA
varComp(Y~trt, dat, ~ibs(SNP))  ## 29 observations compared to 30 in vcf0
dat$SNP[2,1]=NA
varComp(Y~trt, dat, ~ibs(SNP))  ## still, 29 observations, as ibs handles sporadic NA
dat$SNP[3, ]=NA
varComp(Y~trt, dat, ~ibs(SNP))  ## 28 observations, as no genotype for the 3rd obs



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression}
\keyword{ smooth}
