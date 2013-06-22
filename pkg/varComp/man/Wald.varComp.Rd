\name{Wald}
%\alias{KR}
\alias{KR.varComp}
%\alias{satterth}
\alias{satterth.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Wald-type tests for fixed effect parameters
}

\description{
\code{satterth.varComp} computes the denominator degrees of freedom for testing fixed-effect parameters under the variance component model using the Satterthwaite-type method. \code{KR.varComp} computes the denominator degrees of freedom and tests fixed-effect parameters under the variance component model using the Kenward-Roger method. These functions are most conveniently used by \code{\link{fixef.varComp}} and \code{\link{anova.varComp}}. 
}

\usage{
satterth.varComp(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
KR.varComp(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
%\method{satterth}{varComp}()
%\method{KR}{varComp}(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An object of class \code{varComp}.
}
  \item{Lmat}{
	A matrix specifying the linear combinations of fixed effect parameters to be tested for nullity. Each row is a linear combination. 
}
  \item{Vbet}{
	An optional matrix of variance-covariance matrix of fixed-effect parameter estimates.
}
  \item{svd.VLbet}{
	An optional singular value decomposition of \code{Lmat} times \code{Vbet}.
}
  \item{X}{
	Optional fixed-effect design matrix
}
  \item{K}{
	Optional list of \eqn{\mathbf{K}}{K} matrices (the same as in \code{\link{varComp.test}}).
}
  \item{V}{
	Optional covariance matrix of the response variable.
}
  \item{\dots}{
	Place holder, not used. 
}
}

\value{
A numeric scalar of denominator degree of freedom. For \code{satterth.varComp}, the \code{'individual.df'} attribute will contain a vector of numerator degrees of freedom for each row of \code{Lmat}. For \code{KR.varComp}, the result will contain the following attributes: 
\itemize{
\item \code{numDF}:  The numerator degree of freedom, i.e., the rank of \code{Lmat}. 
\item \code{Scale}: A positive numeric scalar, to be multiplied to raw \eqn{F}-statistics before calculating p-values. 
\item \code{F value}: The scaled \eqn{F}-statistic, after adjusting for the variance estimate of fixed effect parameter estimates, and multiplied by the scaling factor. 
\item \code{Pr(>F)}: A numeric scalar of Kenward-Roger \eqn{p}-value. 
\item \code{vcov.beta}: The adjusted variance-covariance estimate of fixed effect parameter estimates. 
}
}
\references{
Michael G. Kenward and James H. Roger (1997) Small sample inference for fixed effects from restricted maximum likelihood. \emph{Biometrics} 53, 983--997

Waseem S. Alnosaier (2007) Kenward-Roger Approximate \eqn{F} Test for Fixed Effects in Mixed Linear Models. \emph{Oregon State University Department of Statistics Ph.D. Thesis}


Fai and Cornelius (1996) Approximate \eqn{F}-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments. \emph{Journal of Statistical Computation and Simulation} 54, 363-378
}
\author{
	Long Qu
}
\seealso{
	\code{\link{fixef.varComp}}, \code{\link{anova.varComp}}
}
\examples{
library(nlme)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)
anova(lmef)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
KR.varComp(vcf, matrix(c(0,1), 1))  # test Source effect
satterth.varComp(vcf, matrix(c(0,1), 1))  # d.f. for testing Source effect

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
