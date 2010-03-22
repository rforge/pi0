\name{print.subt}
\alias{print.subt}
\alias{summary.subt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print the Summary of the Subsampling T-tests }
\description{
This function print out the summary of an \code{subt} object, typically from 
calling the function \code{\link{subt}}.
}
\usage{
\method{print}{subt}(x,...)
\method{summary}{subt}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,object}{ the \code{subt} object, for which to print summaries. }
  \item{...}{ ignored. }
}
\value{
an \code{invisible(NULL)}, used for side effects only. 
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }
\seealso{ \code{\link{subt}} }
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf)
}
data(simulatedSubt)
print(simulatedSubt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ methods }% __ONLY ONE__ keyword per line
