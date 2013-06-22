\name{model.matrix.varComp}
\alias{model.matrix.varComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracting model matrices
}
\description{
This function extracts the fixed-effect design matrix, random-effect design matrix, or the list of variance-covariance matrices whose weighted sum being the variance-covariance matrix of the response variable. 
}
\usage{
\method{model.matrix}{varComp}(object, what = c("fixed", "random", "varcov", "X", "K", "Z"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	A \code{varComp} object
}
  \item{what}{
	A character vector (only the first element will be used) specifying what kind of design matrix is requested. See details. 
}
  \item{\dots}{
	Not used. 
}
}
\details{
\code{"fixed"} and \code{"X"} are equivalent, requesting the fixed effect design matrix.

\code{"random"} and \code{"Z"} are equivalent, requesting the random effect design matrix. Note that this is an equivalent version of the design matrices such that the \code{\link{tcrossprod}} will be the contribution to the marginal variance-covariance. This is not necessarily the one computed directly from the \code{random} argument passed to \code{\link{varComp}}. These are actually computed from \code{\link{cholRoot}} of \eqn{\mathbf{K}}{K} matrices. 

\code{"varcov"} and \code{"K"} are equivalent, requesting the contribution of each random effect to the marginal correlation matrix. These are not necessarily the same value passed to the \code{varcov} argument of \code{\link{varComp}}, because the input value will be treated as the "\eqn{G}" matrices when \code{random} is not missing, but the result here will always be "\eqn{\mathbf{K}}{K}" matrices. See \code{\link{varComp}} for notations. 

}
\value{
If \code{what="fixed"} or \code{"X"}, a single numeric matrix of fixed-effect design matrix. 

Otherwise, a list of requested matrices. 
}
\author{
	Long Qu
}
\note{
See details on possible confusions. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{varComp}}
}
\examples{
library(nlme)
data(Oxide)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
model.matrix(vcf, 'fixed')
model.matrix(vcf, 'random')
model.matrix(vcf, 'varcov')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }