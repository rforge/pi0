\name{p.value.varComp.test}
\alias{p.value}
\alias{p.value.default}
\alias{p.value.htest}
\alias{p.value.varComp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting p-values
}
\description{
\code{p.value} is a generic function for extracting p-values of tests. 
}
\usage{
p.value(x)
\method{p.value}{default}(x)
\method{p.value}{htest}(x)
\method{p.value}{varComp.test}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An R object from which to extract p-values. 
}
}
\details{
By default, \code{p.value} extract the \code{$p.value} component from an \code{htest} object. 

The \code{p.value.varComp.test} function extracts p-values from all tests in the object. 
}
\value{
A numeric vector of p-value(s)
}
\author{
Long Qu
}
\seealso{
\code{\link{varComp.test}}
}
\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
coef(vcf, 'varComp')
p.value(varComp.test(vcf))  ## test both Lot and Wafer%in%Lot
p.value(varComp.test(vcf, null=1))  ## test Wafer%in%Lot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest}
