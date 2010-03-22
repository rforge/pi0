\name{plot.subt}
\alias{plot.subt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 3D Scatter Plot of Subsample Sizes and P-value Density at One. }
\description{
  This function generates a 3d scatter plot for objects of class \code{subt}, 
  typically resulting from calling the function \code{subt}.
}
\usage{
\method{plot}{subt}(x,y,rgl=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the \code{subt} object. }
  \item{y}{ the same as \code{rgl}. If not missing, it overrides \code{rgl}. }
  \item{rgl}{ logical, specifying whether or not the \code{\link[rgl]{plot3d}} in the \pkg{rgl} package is used. 
	If \code{FALSE}, \code{\link[scatterplot3d]{scatterplot3d}} will be used. }
  \item{...}{ other arguments to be passed to either \code{\link[rgl]{plot3d}}, or \code{\link[scatterplot3d]{scatterplot3d}}. }
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
\author{ Long Qu [\email{longor@iastate.edu}] }

\seealso{ \pkg{rgl},\code{\link{subt}} }
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf)
plot(simulatedSubt)
plot(simulatedSubt,FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dynamic }
\keyword{ hplot }% __ONLY ONE__ keyword per line
\keyword{ aplot }% __ONLY ONE__ keyword per line
\keyword{ iplot }% __ONLY ONE__ keyword per line
\keyword{ methods }% __ONLY ONE__ keyword per line
