\name{vcov.varComp}
\alias{vcov.varComp}
\title{Extracting Variance-Covariance Matrices
}
\description{
Extracting (approximate) variance-covariance matrices for fixed-effect parameters, variance components, ratios of variance components to error variance, or the response variable. 
}
\usage{
\method{vcov}{varComp}(object, what = c("fixed", "beta", "random", "varComp", 
	"var.ratio", "tau", "response", "Y"), drop = TRUE, beta.correction=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \code{varComp}.
}
  \item{what}{
	A character vector (only the first element will be used) specifying what variance-covariance matrices are requested. \code{"fixed"} or \code{"beta"} request approximate variance-covariance of fixed-effect parameters (see details). \code{"random"} or \code{"varComp"} request the approximate variance-covariance matrix of variance components computed from the expected information matrix. \code{"var.ratio"} or \code{"tau"} requests approximate variance-covariance matrix of ratio of variance components to the error variance computed from the observed information matrix. \code{"response"} or \code{"Y"} request the marginal variance of the response variable computed from the plug-in estimate. 
}
  \item{drop}{
	A logical scalar, indicating whether zero variance components should be dropped from the results. 
}
	\item{beta.correction}{
	A logical scalar, only applicable when \code{what='beta'}, indicating whether the variance-covariance matrix for fixed effect estimates is corrected according to Kackar and Harville (1984).
}
  \item{\dots}{
	Place holder.
}
}
\details{
For fixed-effect parameters, the results is the plug-in estimate variance of generalized least squares estimates when \code{beta.correction=FALSE}; Otherwise, the Kackar and Harville (1984) correction will be used (default).  For ratios of variance components to error variance, the result is the Hessian matrix. For response variable, the result is the plug-in estimate of the marginal variance. For variance components, the result is the plug-in estimate of inverse expected information matrix from the restricted likelihood.
}
\value{
A numeric matrix of the requested variance-covariance.
}
\references{
Raghu N. Kackar and David A. Harville (1984) Approximations for standard errors of estimators of fixed and random effect in mixed linear models. \emph{Journal of the American Statistical Association}  79, 853--862
}
\author{
 Long Qu
}
\seealso{
\code{\link{varComp}} for the varComp object; 
\code{\link{KR.varComp}} for testing fixed effect parameters accounting for uncertainty in variance parameter estimates. 
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
