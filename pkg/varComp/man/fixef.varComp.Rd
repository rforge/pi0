\name{fixef.varComp}
\alias{fixef}
\alias{fixef.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracting and testing fixed effect parameters
}
\description{
\code{fixef} extracts and tests fixed effect parameters, as well as providing standard errors, etc.
}
\usage{
\method{fixef}{varComp}(object, Lmat, alpha = 0.05, test = "Satterthwaite", ...)
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
  \item{\dots}{
	Not used.
}
}

\value{
A numeric matrix with linear combination of fixed effect parameters, standard errors, degrees of freedom, test statistics, the p-values. Also, the \code{attr( , 'Overall')} attribute is set to the overall test of all requested linear combination of parameters. 
}

\references{
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
data(Orthodont)
vcf = varComp(distance~age*Sex, Orthodont, ~Subject)
fixef(vcf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ regression }

