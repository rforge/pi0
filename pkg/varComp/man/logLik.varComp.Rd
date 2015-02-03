\name{logLik.varComp}
\alias{logLik.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting Profiled Restricted Log Likelihood 
}
\description{
Extracting maximized profiled restricted log likelihood (PREML) from a \code{varComp} object
}
\usage{
\method{logLik}{varComp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	A \code{varComp} object
}
  \item{\dots}{
	Not used
}
}
\details{
	The likelihood value is the profiled restricted log likelihood. The actual value depends on the residual contrast being chosen, but the \emph{location} of maximum does not. Thus it is only comparable when the same residual contrast is used. This value is not to be compared to results reported using functions in other packages, e.g., \code{\link[nlme:lme]{nlme::lme}} or \code{SAS}. The degree of freedom reported is the number of variance components not on the boundary of parameter space, whereas the decision is based on \code{object$control$boundary.eps}. 
}
\value{
	A numeric value of class \code{logLik}.
}
\references{
	Qu L, Guennel T, Marshall SL. (2013) Linear Score Tests for Variance Components in Linear Mixed Models and Applications to Genetic Association Studies. \emph{Biometrics}, Volume 69, Issue 4,  pages 883--892. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link{varComp}}. 
}
\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
logLik(vcf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
