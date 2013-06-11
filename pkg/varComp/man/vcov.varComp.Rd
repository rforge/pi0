\name{vcov.varComp}
\alias{vcov.varComp}
\title{Extracting Variance-Covariance Matrices
}
\description{
Extracting (approximate) variance-covariance matrices for fixed-effect parameters, variance components, ratios of variance components to error variance, or the response variable. 
}
\usage{
\method{vcov}{varComp}(object, what = c("fixed", "beta", "random", "varComp", "var.ratio", "tau", "response", "Y"), drop = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \code{varComp}
}
  \item{what}{
	A character vector (only the first element will be used) specifying what variance-covariance matrices are requested. "fixed" or "beta" request approximate variance-covariance of fixed-effect parameters. "random" or "varComp" request the approximate variance-covariance matrix of variance components. "var.ratio" or "tau" requests approximate variance-covariance matrix of ratio of variance components to the error variance. "response" or "Y" request the marginal variance of the response variable.
}
  \item{drop}{
	A logical scalar, indicating whether zero variance components should be dropped from the results. 
}
  \item{\dots}{
	Place holder.
}
}
\details{
For fixed-effect parameters, the results is the plug-in estimate variance of generalized least squares estimates. For ratios of variance components to error variance, the result is the Hessian matrix. For response variable, the result is the plug-in estimate of the marginal variance. For variance components, the result is the plug-in estimate of inverse information matrix from the restricted likleihood.
}
\value{
A numeric matrix of the requested variance-covariances.
}
\author{
 Long QU
}
\seealso{
\code{\link{varComp}} for the varComp object; 
\code{\link{KR}} for testing fixed effect parameters accounting for uncertainty in variance parameter estimates. 
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
