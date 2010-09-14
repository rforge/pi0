\name{pi0-package}
\alias{pi0-package}
\alias{pi0}
\docType{package}
\title{
Estimating the proportion of true null hypotheses and False Discovery Rates
}
\description{
This package implements method(s) to (approximately unbiasedly) estimate the proportion of true null hypotheses, i.e., the pi0, when a very large number of hypotheses are simultaneously tested, especially for the purpose of (local) false discovery rate control for microarray data. 
}
\details{
\tabular{ll}{
Package: \tab pi0\cr
Type: \tab Package\cr
Version: \tab 1.1-101\cr
Date: \tab 2010-09-14\cr
License: \tab GPL version 2 or newer\cr
}
\itemize{
\item{\code{\link{subt}}}{Subsampling a microarray data set, do t-test for each gene, and estimate p-value density at 1 for each subsample. }
\item{\code{\link{extrp.pi0}}}{Extrapolate the p-value density at 1 over subsample sizes to estimate the proportion of true null hypotheses. }
\item{\code{\link{fdr}}}{Estimate false discovery rate based on p-values and a given estimate of the proportion of true null hypotheses. }
\item{\code{\link{subex}}}{A wrapper that automates \code{\link{subt}},\code{\link{extrp.pi0}},and \code{\link{fdr}}.}
\item{\code{\link{combn2R}}}{ Generating a sample of combinations by choosing m1 out of n1 and m2 out of n2 simultaneously.}
\item{\code{\link{matrix.t.test}}}{Apply a t-test to each row or column of a matrix. }
\item{\code{\link{lastbin}}}{Estimate p-value density at 1 based on a histogram. }
\item{\code{\link{parncpt}}}{Parametrically estimate the distribution of noncentrality parameters. }
\item{\code{\link{nparncpt}}}{Nonparametrically estimate the distribution of noncentrality parameters. }
}
}
\author{
Long Qu \email{longor@iastate.edu}

Maintainer: Long Qu \email{longor@iastate.edu}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\keyword{ package }
\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ nonlinear }
\keyword{ models }
\keyword{ optimize }
\keyword{ smooth }
\keyword{ distribution }
\seealso{
\code{\link{subex}}, \code{\link{subt}}, \code{\link{extrp.pi0}}, \code{\link{fdr}}, \code{\link{combn2R}}
}
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
