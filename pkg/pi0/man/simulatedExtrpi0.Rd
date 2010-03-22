\name{simulatedExtrpi0}
\alias{simulatedExtrpi0}
\docType{data}
\title{ A Simulated 'extrpi0' Object}
\description{
This is the result from calling \code{\link{extrp.pi0}} on \code{\link{simulatedSubt}}.
}
\usage{data(simulatedExtrpi0)}
\format{
an object of class \code{extrpi0}.
}
\details{
This is the result from calling \cr
\code{
data(simulatedSubt)} \cr \code{
simulatedExtrpi0=extrp.pi0(simulatedSubt)
}
}

\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\seealso{\code{\link{extrp.pi0}}}
\examples{
data(simulatedExtrpi0)
print(simulatedExtrpi0)
\dontrun{
plot(simulatedExtrpi0)
}
}
\keyword{datasets}
