\name{simulatedSubex}
\alias{simulatedSubex}
\docType{data}
\title{ A Simulated 'subex' Object }
\description{
This is the result from calling \code{\link{subex}} on \code{\link{simulatedDat}}.
}
\usage{data(simulatedSubex)}
\format{
An object of class \code{\link{subex}}
}
\details{
This is the result from calling \cr
\code{
data(simulatedDat)}\cr \code{
simulatedSubex=subex(simulatedDat,balanced=FALSE,max.reps=Inf,plotit=FALSE)
}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\seealso{\code{\link{subex}}}
\examples{
data(simulatedSubex)
print(simulatedSubex)
\dontrun{
plot(simulatedSubex)
}
}

\keyword{datasets}
