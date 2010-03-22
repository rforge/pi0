\name{plot.subex}
\alias{plot.subex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the P-value, Q-values, and the Regression Surface }
\description{
This function plots the p-value, q-values, and the regression surface for an object of class \code{\link{subex}},
typically from calling the function \code{\link{subex}}. 
}
\usage{
\method{plot}{subex}(x,y,rgl = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the \code{subex} object }
  \item{y}{ the same as \code{rgl}. If not missing, it overrides \code{rgl}. }
  \item{rgl}{ logical, specifying whether or not the \pkg{rgl} package is used for making better 3D interactive graphs. }
  \item{...}{ other arguments to be passed to \code{\link[rgl]{persp3d}}, or \code{\link{persp}} if \pkg{rgl} package is not available. }
}
\details{
	Two plots will be generated. The first one is a histogram of p-values in the blue color. 
	A horizontal blue line is added indicating the height of \eqn{\pi_0}{pi0}, i.e., the proportion
	of true null hypotheses. This histogram is overlaid with a red line of FDRs, indicating the corresponding 
	q-value for each p-value. The right most end, i.e., the q-value corresponding to a p-value of 1, 
	is the also of height \eqn{\pi_0}{pi0}.\cr
	The other plot is the same as the \code{\link{plot.extrpi0}}.
}
\value{
an \code{invisible(NULL)}, used for side effects only.
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }
\note{
  Because this function will call \code{\link{plot.extrpi0}}, a warning will be generated when \pkg{rgl} 
  package is not available. 
}
\seealso{ \code{\link{plot.extrpi0}}, \pkg{rgl} ,\code{\link{subex}} }
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubex' object in this package generated
simulatedSubex=subex(simulatedDat,balanced=FALSE,max.reps=Inf,plotit=FALSE)
plot(simulatedSubex)
plot(simulatedSubex,FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dynamic }
\keyword{ hplot }% __ONLY ONE__ keyword per line
\keyword{ aplot }% __ONLY ONE__ keyword per line
\keyword{ iplot }% __ONLY ONE__ keyword per line
\keyword{ methods }% __ONLY ONE__ keyword per line
