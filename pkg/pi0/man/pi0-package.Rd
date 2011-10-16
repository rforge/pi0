\name{pi0-package}
\alias{pi0-package}
\alias{pi0}
\docType{package}
\title{
Estimating the proportion of true null hypotheses and False Discovery Rates
}
\description{
This package implements method(s) to (approximately unbiasedly) estimate the proportion of true null hypotheses, i.e., the pi0, when a very large number of hypotheses are simultaneously tested, especially for the purpose of (local) false discovery rate control for microarray data. It also contains functions to estimate the distribution of noncentrality parameters from a large number of parametric tests. 
}
\details{
\tabular{ll}{
Package: \tab pi0\cr
Type: \tab Package\cr
Version: \tab 1.3-190\cr
Date: \tab 2011-10-15\cr
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
\item{\code{\link{sparncpt}}}{Semiparametrically estimate the distribution of noncentrality parameters.} 
\item{\code{\link{nparncpp}}}{Nonparametric estimate of the distribution of absolute noncentrality parameters from a large number of p-values. }
\item{\code{\link{CBUM}}}{(Censored) Beta-Uniform mixture model for p-values. }
\item{\code{\link{znormix}}}{Normal mixture model for z-scores. }
}
}
\author{
Long Qu \email{lqu@biostatsolutions.com}

Maintainer: Long Qu \email{lqu@biostatsolutions.com}
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.

Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)

Ruppert, Nettleton, Hwang. 2007. Exploring the Information in $p$-values for the Analysis and Planning of Multiple-test Experiments. Biometrics. 63. 483-495.

G.J. McLachlan, R.W. Bean and L. Ben-Tovim Jones. (2006) A Simple implementation of a normal mixture approach to differential gene expression in multiclass microarrays. Bioinformatics, 22(13):1608-1615. 

Anastasios Markitsis and Yinglei Lai (2010) A censored beta mixture model for the estimation of the proportion of non-differentially expressed genes. Bioinformatics 26(5):640-646.

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
\code{\link{subex}}, \code{\link{subt}}, \code{\link{extrp.pi0}}, \code{\link{fdr}}, \code{\link{combn2R}}, 
\code{\link{nparncpt}}, \code{\link{parncpt}}, \code{\link{sparncpt}}, \code{\link{nparncpp}}
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
data(simulatedSubex); print(simulatedSubex)
## parametric, nonparametric, semiparametric estimate of 
## noncentrality parameter distribution from t-statistics
data(simulatedTstat)
(npfit=nparncpt(tstat=simulatedTstat, df=8, plotit=FALSE)); 
(pfit=parncpt(tstat=simulatedTstat, df=8, zeromean=FALSE)); 
(pfit0=parncpt(tstat=simulatedTstat, df=8, zeromean=TRUE)); 
(spfit=sparncpt(npfit,pfit));
}
