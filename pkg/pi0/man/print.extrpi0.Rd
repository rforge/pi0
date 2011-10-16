\name{print.extrpi0}
\alias{print.extrpi0}
\alias{summary.extrpi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print the Summary of the Extrapolation }
\description{
This function print out the summary of an \code{extrpi0} object, typically from 
calling the function \code{\link{extrp.pi0}}.
}
\usage{
\method{print}{extrpi0}(x,...)
\method{summary}{extrpi0}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,object}{ the \code{extrpi0} object, for which to print summaries. }
  \item{...}{ ignored. }
}
\value{
an \code{invisible(NULL)}, used for side effects only. 
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{lqu@biostatsolutions.com} }
\seealso{ \code{\link{extrp.pi0}} }
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf)
## this is how the 'simulatedExtrpi0' data set in this package generated
simulatedExtrpi0=extrp.pi0(simulatedSubt)
}
data(simulatedExtrpi0)
print(simulatedExtrpi0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ methods }% __ONLY ONE__ keyword per line
