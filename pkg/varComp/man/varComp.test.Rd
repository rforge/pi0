\name{varComp.test}
\alias{varComp.test}
\alias{varComp.test.formula}
\alias{varComp.test.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests for Nullity of Variance Components
}
\description{
These are the user-interface functions to test the nullity of variance components, using linear score tests, projected quadratic score tests, or restricted likelihood ratio tests. 
}
\usage{
varComp.test(object, ...)
\method{varComp.test}{formula}(object, data, random1, varcov1, 
random2, varcov2, fit.control, test.control, ...)
\method{varComp.test}{varComp}(object, object2, additional.varcov, null, 
test = "LinScore", control = varCompTest.control(test), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, object2}{
	An R object. In \code{varComp.test.formula}, this is the fixed-effect formula. In \code{varComp.test.varComp} these are the \code{\link{varComp}} objects specifying the null model and the alternative models. \code{object2} could be missing, but then either \code{additional.varcov} or \code{null} needs to be given.
}
  \item{data}{
	The same as in \code{\link{varComp}}.
}
  \item{random1, random2}{
	The random effect formulas for the first and the second models, respectively. Specification is the same as in \code{\link{varComp}}. 
}
  \item{varcov1, varcov2}{
	The lists of variance-covariance matrices for the first and the second models, respectively. Specification is the same as in \code{\link{varComp}}. 
}
  \item{fit.control}{
	A list of arguments of class \code{\link{varComp.control}}. 
}
  \item{test.control, control}{
	A list of arguments of class \code{\link{varCompTest.control}}. 
}
  \item{additional.varcov}{
	An optional list of variance-covariance matrices. It is only used when \code{object2} is missing. When this is given, \code{object} will be treated as the null model. And the test of interest is the nullity of variance components corresponding to these \code{additional.varcov}.  \code{null} needs to be missing when \code{additional.varcov} is given. 
}
  \item{null}{
	An optional integer vector, indicating the variance components in the null model. It is only used when \code{object2} is missing, and is incompatible with \code{additional.varcov}. When \code{null} is given, \code{object} is treated as the alternative model. And the test of interest is the nullity of variance components \emph{not} included in \code{null}. For example, if \code{object} has two variance components (other than the error variance). Setting \code{null=1L} will test the nullity of the second components, whereas setting \code{null=integer(0L)} (the default) will test the simultaneous nullity of both components. 
}
  \item{test}{
	A character vector, specifying the tests to be performed. This does not need to be a vector of length 1. Available choices are given in \code{\link{varCompTests}}. Currently, this includes \code{"LinScore"}, \code{"VM03"}, \code{"SS95"}, \code{"HP01"}, and \code{"RLRT"}. However, not all test methods are applicable to all models/hypotheses. 
}
  \item{\dots}{
	Place holder.
}
}
\details{
The formula interface \code{varComp.test.formula} is primarily used to test two variance component models with the same fixed-effect formula. The two models need to be nested. 

In the \code{varComp.test.varComp} interface, exactly \emph{one} of \code{object2}, \code{additional.varcov} and \code{null} needs to be given. If \code{object2} is given, two-model comparison is performed against \code{object}. If \code{additional.varcov} is given, \code{object} is treated as the null and \code{additional.varcov} is treated as the additional variance components to be tested.  If \code{null} is given, \code{object} is treated as the alternative and the components not in \code{null} will be tested. 
}
\value{
A list of class \code{varComp.test} of the same length as \code{test}. Each component is either an object of \code{htest} or a list of \code{htest} objects, depending upon whether multiple methods are requested under the same \code{test}.  See \code{\link{varCompTest.control}} for more details on controlling the tests to be performed. 
}
\references{

	Qu L, Guennel T, Marshall SL. (2013) Linear Score Tests for Variance Components in Linear Mixed Models and Applications to Genetic Association Studies. \emph{Biometrics}, Volume 69, Issue 4,  pages 883--892. 
	
	Daniel B. Hall and Jens T. Praestgaard (2001) Order-restricted tests for homogeneity in generalised linear and nonlinear mixed models. \emph{Biometrika}, 88, 739--751. 
	
	Mervyn J. Silvapulle and Paramsothy Silvapulle (1995) A score test against one-sided alternatives. \emph{Journal of the American Statistical Association}, 90, 342--349.
	
	Geert Verbeke and Geert Molenberghs (2003) The use of score tests for inference on variance components. \emph{Biometrics}, 59, 254--262. 
	
	Ciprian M. Crainiceanu and David Ruppert (2003) Likelihood ratio tests in linear mixed models with one variance component. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 66, 165--185.
	
	Greven S, Crainiceanu CM Kuchenhoff H and Peters A (2008) Restricted Likelihood Ratio Testing for Zero Variance Components in Linear Mixed Models. \emph{Journal of Computational and Graphical Statistics}, 17, 870--891.
}
\author{
Long Qu
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{varCompTest.control}}, \code{\link{varComp}}
}
\examples{
### Oxide/Semiconductor data example
library(nlme)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
VarCorr(lmef)
coef(vcf, 'varComp') ## same values as above
varComp.test(vcf)  ## test against linear model
varComp.test(vcf, null=1)  ## test against model with only Lot random effect


### Genetics example	
trt=gl(2, 15)
set.seed(2340)
dat=data.frame(trt=trt)
dat$SNP=matrix(sample(0:2, 120, replace=TRUE), 30)
dat$Y  = as.numeric(trt)+rnorm(30) + dat$SNP\%*\%rnorm(4)
(vcf0 = varComp(Y~trt, dat, ~ibs(SNP)))
(vcf1 = varComp(Y~trt, dat, ~ibs(SNP):trt)) ## two variance components

varComp.test(vcf0)
varComp.test(vcf1)
varComp.test(vcf0, vcf1)
varComp.test(vcf1, null=1)
varComp.test(Y~trt, dat, random1=~ibs(SNP), random2=~trt:ibs(SNP))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
