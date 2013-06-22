\name{coef.varComp}
\alias{coef.varComp}
%\alias{fixef}
%\alias{fixed.effects}
%\alias{fixef.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Obtaining parameter estimates from a varComp object
}
\description{
\code{coef} computes and/or extracts fixed-effect parameters, variance-component parameters, or ratios of variance component parameters to the error variance. %\code{fixef} and \code{fixed.effects} are equivalent to \code{coef} with \code{what='fixed'}.
}
\usage{
\method{coef}{varComp}(object, what = c("fixed", "beta", "random", "varComp", "var.ratio", "tau"), ...)
%\method{fixef}{varComp}(object, ...)
%\method{fixed.effects}{varComp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \code{varComp}
}
  \item{what}{
	Character vector (only the first component will be used) specifying what parameters are requested. See details. 
}
  \item{\dots}{	Not used.
}
}
\details{
\code{"fixed"} (default) and \code{"beta"} are equivalent, requesting fixed-effect parameters to be returned. 

\code{"random"} and \code{"varComp"} are equivalent, requesting variance components (including the error variance) to be returned. 

\code{"var.ratio"} and \code{"tau"} are equivalent, requesting the ratio of variance components to the error variance to be returned. 
}
\value{
	A named numeric vector of requested parameter estimates.
}

\author{
Long Qu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{varComp}}
}
\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
coef(vcf, 'varComp') 
coef(vcf )  ## same as coef(vcf, 'fixed')
coef(vcf, 'var.ratio')

%fixef(vcf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ models }
