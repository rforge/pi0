\name{Wald}
\alias{KR}
\alias{KR.varComp}
\alias{satterth}
\alias{satterth.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Wald-type tests for fixed effect parameters
}

\description{
\code{satterth} tests fixed-effect parameters under the variance component model using the Satterthwaite-type method. \code{KR.varComp} tests fixed-effect parameters under the variance component model using the Kenward-Roger method
}

\usage{
satterth(object, ...)
KR(object, ...)
\method{satterth}{varComp}(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
\method{KR}{varComp}(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An object of class \code{varComp}
}
  \item{Lmat}{
	A matrix specifying the linear combinations of fixed effect parameters to be tested for nullity. Each row is a linear combination. 
}
  \item{Vbet}{
	An optional matrix of variance-covariance matrix of fixed-effect parameter estimates
}
  \item{svd.VLbet}{
	An optional singular value decomposition of \code{Lmat} times \code{Vbet}.
}
  \item{X}{
	Optional fixed-effect design matrix
}
  \item{K}{
	Optional list of K matrices (the same as in \code{\link{varComp.test}}.
}
  \item{V}{
	Optional covariance matrix of the response variable.
}
  \item{\dots}{
	Place holder, not used. 
}
}

\value{
A numeric scalar of p-value
}
\references{
Michael G. Kenward and James H. Roger (1997) Small sample inference for fixed effects from restricted maximum likelihood. \emph{Biometrics} 53, 983--997

Waseem S. Alnosaier (2007) Kenward-Roger Approximate F Test for Fixed Effects in Mixed Linear Models. \emph{Oregon State University Department of Statistics Ph.D. Thesis}


Fai and Cornelius (1996) Approximate F-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments. \emph{Journal of Statistical Computation and Simulation} 54, 363-378
}
\author{
	Long Qu
}
\seealso{
	\code{\link{satterth}}
}
\examples{
library(nlme)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)
anova(lmef)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
KR(vcf, matrix(c(0,1), 1))  # test Source effect
satterth(vcf, matrix(c(0,1), 1))  # test Source effect

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
