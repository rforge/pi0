\name{anova.varComp}
\alias{anova.varComp}
\alias{KR.varComp}
\alias{satterth.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ANOVA-type analysis of fixed effect parameters
}
\description{
\code{anova.varComp} tests for fixed effect contrasts, as well as providing standard errors, etc. \code{satterth.varComp} computes the Satterthwaite-type denominaor degrees of freedom. \code{KR.varComp} computes the Kenward-Roger denominaotr degrees of freedom, adjusted stanard errors, and scale factors applied to the F-statistic. 
}
\usage{
\method{anova}{varComp}(object, Lmat, alpha = 0.05, test = "KR", ...)
satterth.varComp(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
KR.varComp(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An object of class \code{varComp}
}
  \item{Lmat}{
	A matrix specifying the linear combinations of fixed effect parameters to be tested for nullity. Each row is a linear combination. 
}
  \item{alpha}{
	The requested level of test. 
}
  \item{test}{
	A character scalar requesting the type of the tests to be performed. Currently, only \code{"Satterthwaite"} is supported. 
}
  \item{Vbet}{
	An optional matrix of variance-covariance matrix of fixed-effect parameter estimates \code{without} any correction. 
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
	Not used.
}
}

\value{
A numeric matrix with linear combination of fixed effect parameters, standard errors, degrees of freedom, test statistics, the p-values. Also, the \code{attr( , 'Overall')} attribute is set to the overall test of all requested linear combination of parameters. 

\code{satterth.varComp} returns a named numeric scalar, which is the denominaor degrees of freedom. It also has the \code{individual.df} attribute, which is a vector of denominator degrees of freedom for each row of \code{Lmat}. 

\code{KR.varComp} returns a named numeric scalar, which is the denominaor degrees of freedom. It also has the following attributes: 
\itemize{
	\item \code{numerator.df}: A numeric scalar of the row rank of \code{Lmat}, i.e., the numerator degrees of freedom. 
	\item \code{scale}: A numeric scalar of the scale factor to be multiplied to the raw F-statistic. 
	\item \code{raw.F}; A numeric scalar of the raw F-statistic. 
	\item \code{p.value}: A numeric scalar of the p-value. 
	\item \code{vcov.beta}: A numeric matrix of the corrected variance-covariance matrix of fixed effect estimates. 
}
}

\references{

Michael G. Kenward and James H. Roger (1997) Small sample inference for fixed effects from restricted maximum likelihood. \emph{Biometrics} 53, 983--997

Waseem S. Alnosaier (2007) Kenward-Roger Approximate F Test for Fixed Effects in Mixed Linear Models. \emph{Oregon State University Department of Statistics Ph.D. Thesis}


Fai and Cornelius (1996) Approximate F-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments. \emph{Journal of Statistical Computation and Simulation} 54, 363-378
}

\author{
Long Qu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{satterth}}
}

\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
fixef(vcf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ regression }

