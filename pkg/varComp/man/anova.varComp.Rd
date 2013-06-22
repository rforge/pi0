\name{anova.varComp}
\alias{anova.varComp}
\alias{fixef}
\alias{fixed.effects}
\alias{fixef.varComp}
\alias{print.varCompFixEf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ANOVA-type analysis of fixed effect parameters
}
\description{
\code{fixef.varComp} and \code{anova.varComp} test for fixed effect contrasts, as well as providing standard errors, etc. 
}
\usage{
\method{anova}{varComp}(object, ..., test = c("KR", 'Satterthwaite'), L)
\method{fixef}{varComp}(object, Lmat, alpha=.05, test=c('KR', 'Satterthwaite'), ...)
\method{print}{varCompFixEf}(x, ...)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{object,x}{
	An object of class \code{varComp}.
}
  \item{Lmat,L}{
	A matrix specifying the linear combinations of fixed effect parameters to be tested for nullity. Each row is a linear combination. If the column names are non-null, the columns will be re-ordered according to the column names of \code{model.matrix(object)}. In case of a single row matrix, it may also be given as a vector. 
}
  \item{alpha}{
	The requested level of test. 
}
  \item{test}{
	A character scalar requesting the type of the tests to be performed. 
}
  \item{\dots}{
	For \code{anova.varComp}, this may be optional additional \code{varComp} objects. For \code{fixef.varComp} and \code{print.varCompFixEf}, this is not used.
}
}

\value{
For \code{anova.varComp}, this is a data frame with the following columns: 
\itemize{
	\item \code{F value}: The scaled \eqn{F}-statistic; 
	\item \code{Scale}: The multiplicative factor in the scaled \eqn{F}-statistic;
	\item \code{numDF}: The numerator degrees of freedom;
	\item \code{denDF}: The denominator degrees of freedom;
	\item \code{Pr(>F)}: The \eqn{p}-value.
}
When \code{\dots} is missing but \code{L} is given, each row of \code{L} matrix will have a corresponding row in the result, testing the nullity of linear contrasts specified by that row of \code{L}. Additionally, an "Overall" row will be the last row, testing the nullity of linear contrasts of all rows of \code{L}. 

When \code{\dots} is not missing, \code{L} will be ignored and the result will have the same number of rows as the number of \code{varComp} objects given (including the first \code{object} argument). The order of rows is such that the rank of the fixed effect design matrix is non-decreasing, i.e., smaller models first. The results will be the comparison between the current model vs. the model in the preceding row, except for the first row. For the first row, the model being compared is a null model. If the fixed effect design matrix for the model in first row includes the intercept in its column space, the null model includes only an intercept. If the intercept is not in the column space of the first model, the null model will contain no fixed effects. 

When both \code{\dots} and \code{L} are missing, \code{anova.varComp} will test each fixed effect parameters together with a last "Overall" row. However, this last row will remove the intercept from the model, if there is one. See example below. 

For \code{fixef.varComp}, the result is an object with class \code{varCompFixEf}. This is a named numeric vector, providing estimates of  fixed effect contrasts specified by \code{Lmat}. The name will be the row names of \code{Lmat}. It has the \code{anova} attribute, a numeric matrix with the following columns: 
\itemize{
\item \code{Std. Error}: the standard error of estimate. If \code{test='KR'}, this is the adjusted standard error. If \code{test='Satterthwaite'}, this is the plug-in estimate. 
\item \code{Lower}: the lower limit of the \code{1-alpha} confidence interval. 
\item \code{Upper}: the upper limit of the \code{1-alpha} confidence interval. 
\item \code{t value}: the scaled \eqn{t}-statistic testing nullity of this contrast. 
\item \code{Scale}: a multiplicative factor in the scaled \eqn{t}-statistic. 
\item \code{Df}: the degrees of freedom for the \eqn{t}-statistic. 
\item \code{Pr(>|t|)}: the \eqn{p}-value.
}
Additionally, the \code{anova} attribute also has the \code{Overall} attribute, which is a single-row matrix with the same columns as the results of \code{anova.varComp}, testing for the overall nullity of linear contrasts specified by \code{Lmat}.

\code{print.varCompFixEf} prints the result from \code{fixef.varComp} and return the argument invisibly.
}

\references{

Michael G. Kenward and James H. Roger (1997) Small sample inference for fixed effects from restricted maximum likelihood. \emph{Biometrics} 53, 983--997

Waseem S. Alnosaier (2007) Kenward-Roger Approximate \eqn{F} Test for Fixed Effects in Mixed Linear Models. \emph{Oregon State University Department of Statistics Ph.D. Thesis}


Fai and Cornelius (1996) Approximate \eqn{F}-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments. \emph{Journal of Statistical Computation and Simulation} 54, 363-378
}

\author{
Long Qu
}
\section{Warning }{
%% ~Make other sections like Warning with \section{Warning }{....} ~
\enumerate{
\item When neither \code{\dots} nor \code{L} is given for \code{anova.varComp}, the current implementation will test each fixed effect parameter separately. This behavior might be changed in future releases. 
\item When \code{L} is not explicitly given, the current implementation will rely on \code{getOption("digits")} when trying to find the \code{L} matrix that test the difference between adjacent models.  
}
}

\seealso{
\code{\link{satterth.varComp}}, \code{\link{KR.varComp}}
}

\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
vcf0 = varComp(Thickness~1, Oxide, ~Lot/Wafer)
fixef(vcf)
anova(vcf, vcf0)
anova(vcf)
anova(vcf, L=diag(1,2)) ## notice the difference in the last row
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ regression }
\keyword{ print}

