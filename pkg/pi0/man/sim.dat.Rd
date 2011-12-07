\name{sim.dat}
\alias{sim.dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulating a Microarray Data Set }
\description{
This function simulates a two-group comparison microarray data set according to a hierarchical model, 
where the standardized effect sizes across all genes are assumed to be independently and identically 
distributed. This distribution is a two-component mixture. It has probability \eqn{\pi_0}{pi0} of being 
zero; and probability \eqn{1-\pi_0}{1-pi0} of being from another distribution. The observed values 
are simulated independently conditional on the standardized effect sizes.
}
\usage{
sim.dat(G = 10000, pi0 = 0.75, gamma2 = 1, n1 = 5, n2 = n1, 
        errdist = rnorm, effdist = function(g, gamma2) 
        rnorm(g, , sqrt(gamma2)), ErrArgs, EffArgs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{ a numeric positive integer, the number of genes}
  \item{pi0}{ a numeric value between 0 and 1, the proportion of non-differentially expressed genes. }
  \item{gamma2}{ a positive value, which is always the second argument passed to \code{effdist}. 
	If the nonzero standardized effect sizes have a zero normal distribution, this is the variance 
	of this distribution. The larger it is, the larger the mean absolute effects are. }
  \item{n1}{ a positive integer, the sample size in treatment group 1. }
  \item{n2}{ a positive integer, the sample size in treatment group 2. }
  \item{errdist}{ a function, which simulate \code{K} random errors, where \code{K} is the first argument of \code{errdist}. The second argument is always \code{ErrArgs}, if it is not missing. }
  \item{effdist}{ a function, which simulate \code{G1} standardized effect sizes, where \code{G1} is the first argument of \code{effdist}. The second argument is always \code{gamma2}. The third argument is always \code{EffArgs}, if it is not missing. }
  \item{ErrArgs}{a list of additional arguments used by \code{errdist}.}
  \item{EffArgs}{a list of additional arguments used by \code{effdist}.}
}
\details{
The funciton simulates \eqn{G*N}{G*N} errors according to \code{errdist}, where \eqn{N=n_1+n_2}{N=n1+n2}. The results
are organized into a G-by-N matrix. The \eqn{G_1}{G1} standarized effect sizes are simulated according to 
\code{effdist}, controlled by the parameter \code{gamma2}, where \eqn{G_1=round(G* pi_0 }{\code{G1=round(G*pi0)}}. 
Then, each column of the upper-left  G1-by-n1 submatrix were added by the simulated effect sizes. 
}
\value{
a \code{G}-by-\code{(n1+n2)} matrix. 
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu }
\examples{
set.seed(54457704)
## an unusually small data set of 20 genes and 3 samples in each of the two treatment groups. 
dat=sim.dat(G=20, n1=3,n2=3)

set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
