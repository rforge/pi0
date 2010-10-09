\name{MRPP-package}
\Rdversion{1.1}
\alias{MRPP-package}
\alias{MRPP}
\docType{package}
\title{
Multiresponse permutation procedure and its variable importance and selection. 
}
\description{
Implements multiresponse permutation procedure (MRPP), variable importance measure for MRPP, and backward variable selection procedures for the MRPP, applicable to high dimensional data. 
}
\details{
\tabular{ll}{
Package: \tab MRPP\cr
Type: \tab Package\cr
Version: \tab 1.0-138\cr
Date: \tab 2010-10-09\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
\itemize{
\item{\code{\link{mrpp.test}}} {The MRPP, a permutation test applicable to high dimensional data.}
\item{\code{\link{get.dp.dw.kde}}} {Variable importance measure: approximate derivative of p-values to weights.}
\item{\code{\link{get.p.dd.dw}}} {Variable importance measure: permutation p-value using the derivative of Euclidean distance with respect to weights as the distance measure.}
\item{\code{\link{back.search}}} {Backward variable selection for the MRPP.}
\item{\code{\link{bsmrpp}}} {Backward variable selection nested within outer permutations.}
}
}
\author{
Long Qu

Maintainer: Long Qu <longor@iastate.edu>
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. RELATIVE VARIABLE IMPORTANCE AND BACKWARD VARIABLE SELECTION FOR THE MULTIRESPONSE
PERMUTATION PROCEDURE, WITH APPLICATIONS TO HIGH DIMENSIONAL GENOMIC DATA. (submitted)
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%~~ simple examples of the most important functions ~~
}
