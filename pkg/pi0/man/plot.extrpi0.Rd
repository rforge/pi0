\name{plot.extrpi0}
\alias{plot.extrpi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the Estimated Regression Surface }
\description{
  This function plots the regression surface, overlaid with 3D scatter plot, for objects of class \code{extrpi0}, 
  typically resulting from calling the function \code{extrp.pi0}.
}
\usage{
\method{plot}{extrpi0}(x,y,rgl=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the \code{extrpi0} object. }
  \item{y}{ the same as \code{rgl}. If not missing, it overrides \code{rgl}. }
  \item{rgl}{ logical, specifying whether or not the \pkg{rgl} package is used for making better 3D interactive graphs. }
  \item{...}{ other arguments to be passed to either \code{\link[rgl]{persp3d}}, or \code{\link{persp}} if \pkg{rgl} is not available. }
}

\value{
an \code{invisible(NULL)}. Used for side effects only.
}
\note{
	When \pkg{rgl} is not available, a warnging is always generated.
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }

\seealso{ \pkg{rgl},\code{\link{extrp.pi0}} }
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)	
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf) 
## this is how the 'simulatedExtrpi0' data set in this package generated
simulatedExtrpi0=extrp.pi0(simulatedSubt,plotit=FALSE) 
plot(simulatedExtrpi0)
plot(simulatedExtrpi0,FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dynamic }
\keyword{ hplot }% __ONLY ONE__ keyword per line
\keyword{ aplot }% __ONLY ONE__ keyword per line
\keyword{ iplot }% __ONLY ONE__ keyword per line
\keyword{ methods }% __ONLY ONE__ keyword per line
