\name{subex}
\alias{subex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Subsampling-Extrapolation Based Estimation of Proportion of True Null Hypotheses and False Discovery Rates for Microarray Data }
\description{
This function is a wrapper of \code{\link{subt}}, \code{\link{extrp.pi0}} and \code{\link{fdr}}, and is
a ready to use directly on a matrix of microarray data. 
}
\usage{
subex(dat, n1 = round(ncol(dat)/2), n2 = ncol(dat) - n1, 
  f1method = c("lastbin", "qvalue"), 
  max.reps = 20, balanced = FALSE, nparm = c(2, 4), 
  extrpFUN = c("constrOptim", "genoud"), 
  starts = c(pi0 = 0.75, gam2 = 1, a = 0.5, c = 0.5), plotit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a numeric matrix, which is the microarray data. Each row represent a gene, and each column 
                represent a subject. The first \code{n1} columns correspond to the first treatment group; 
                and the rest \code{n2} columns correspond ot the second treatment group. }
  \item{n1}{ a positive integer, the sample size in treatment group 1. }
  \item{n2}{ a positive integer, the sample size in treatment group 2. }
  \item{f1method}{ character, the name of the function used to estimate the p-value density at one. See \code{\link{subt}} for details. }
  \item{max.reps}{ a positive integer, the maximum number of subsamples "per subsample size configuration". See \code{\link{subt}} for details. }
  \item{balanced}{ logical, indicating if only balanced subsamples are generated.  See \code{\link{subt}} for details. }
  \item{nparm}{ either 2 or 4, indicating the number of parameters used in extrapolation.  See \code{\link{extrp.pi0}} for details.  }
  \item{extrpFUN}{ character, specifying the name of the optimization function for nonlinear regression.  See the \code{FUN} argument of \code{\link{extrp.pi0}} for details.  }
  \item{starts}{ a numeric vector of length \code{nparm}, specifying the starting values of optimization.  See \code{\link{extrp.pi0}} for details.  }
  \item{plotit}{ logical, indicating if the extrapolation plot will be produced.  See \code{\link{extrp.pi0}} for details.  }
}
\details{
This function calls \code{\link{subt}}, \code{\link{extrp.pi0}}, \code{\link{matrix.t.test}} and \code{\link{fdr}} 
sequentially to estimate the proportion of true null hypotheses \eqn{\pi_0}{pi0} as well as the false discovery
rates (FDR) based on the estimated \eqn{\pi_0}{pi0}.
}
\value{
an object of class \code{subex}, which is a list 4 components:
\item{pi0}{ a numeric value, giving the estimated \eqn{\pi_0}{pi0}}.
\item{extrp.fit}{ an object of class \code{extrpi0}, the results from calling \code{\link{extrp.pi0}}.}
\item{pvalues}{a numeric vector of length the same as \code{nrow(dat)}, the p-values for each gene. }
\item{qvalues}{a numeric vector of length the same as \code{nrow(dat)}, the q-values for each gene. }


}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }
\note{ 
Plotting using package \pkg{rgl} will be tried. If not available, a warning will be generated. See \code{\link{plot.extrpi0}} for details.
}
\seealso{\code{\link{subt}}, \code{\link{extrp.pi0}}, \code{\link{matrix.t.test}},\code{\link{fdr}},
\code{\link{plot.subex}},\code{\link{print.subex}} }

\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)
## this is how the 'simulatedSubex' object in this package generated
simulatedSubex=subex(simulatedDat,balanced=FALSE,max.reps=Inf,plotit=FALSE) 
plot(simulatedSubex)
}
data(simulatedSubex)
print(simulatedSubex)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ htest }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{iteration}
\keyword{optimize}
