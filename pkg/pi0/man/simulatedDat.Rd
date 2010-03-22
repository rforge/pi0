\name{simulatedDat}
\alias{simulatedDat}
\docType{data}
\title{ A Simulated Microarray Data Set}
\description{
This is the result from calling \code{\link{sim.dat}}
}
\usage{data(simulatedDat)}
\format{
a matrix with 5000 rows and 10 columns. The first 5 columns correspond to treatment 
one, and the rest 5 columns corresond to treatment 2. Each row corresponds to a gene.
}
\details{
This is the result from calling \cr
\code{
set.seed(9992722)}\cr \code{
simulatedDat=sim.dat(G=5000)
}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\seealso{\code{\link{sim.dat}}}
\examples{
data(simulatedDat)
}
\keyword{datasets}
