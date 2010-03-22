\name{simulatedSubt}
\alias{simulatedSubt}
\docType{data}
\title{ A Simulated 'subt' Object}
\description{
This is the result from calling \code{\link{subt}} on \code{\link{simulatedDat}}
}
\usage{data(simulatedSubt)}
\format{
an object of class \code{\link{subt}} with \code{balanced=FALSE} with \code{max.reps=Inf} attributes.
}
\details{
This is the result from calling \cr
\code{
data(simulatedDat)}\cr \code{
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf)
}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\seealso{\code{\link{subt}}}
\examples{
data(simulatedSubt)
print(simulatedSubt)
\dontrun{
plot(simulatedSubt)
}
}
\keyword{datasets}
