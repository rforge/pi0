\name{simulatedTstat}
\alias{simulatedTstat}
\alias{simulatedPval}
\docType{data}
\title{ t-Statistics and p-Values from a Simulated Microarray Data Set}
\description{
These are the results from applying t-tests to the simulated data set \code{\link{simulatedDat}}.
}
\usage{
data(simulatedTstat)
data(simulatedPval)
}
\format{
A numeric vector of two-sample t-statistics of length 5000.
}
\details{
This is the result from calling \cr
\code{
data(simulatedDat)
simulatedTstat=matrix.t.test(simulatedDat,tOnly=TRUE)
simulatedPval=matrix.t.test(simulatedDat)
}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\seealso{
	\code{\link{sim.dat}}, \code{\link{simulatedDat}}
	}
\examples{
data(simulatedTstat)
data(simulatedPval)
}
\keyword{datasets}
