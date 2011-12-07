\name{print.subex}
\alias{print.subex}
\alias{summary.subex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Printing a Summary of Subsampling-Extrapolation Results }
\description{
This function prints out a summary of a \code{\link{subt}} object. 
}
\usage{
\method{print}{subex}(x,...)
\method{summary}{subex}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,object}{ the \code{subex} object, for which to print summaries. }
  \item{...}{ ignored. }
}
\details{
	This function will first print a summary of the corresponding \code{extrpi0} object. 
	Then several quantiles of the p-values and q-values are printed.
}
\value{
an \code{invisible(NULL)}, used only for its side effects. 
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu }
\seealso{ \code{\link{subex}}, \code{\link{extrp.pi0}}, \code{\link{print.extrpi0}}, \code{\link{fdr}}}
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubex' data set in this package generated
simulatedSubex=subex(simulatedDat,balanced=FALSE,max.reps=Inf)
}
data(simulatedSubex)
summary(simulatedSubex)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ methods }% __ONLY ONE__ keyword per line
