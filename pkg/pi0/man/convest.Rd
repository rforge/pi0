\name{convest}
\alias{convest}
\title{Estimate Proportion of True Null Hypotheses}
\description{
Returns an estimate of the proportion of true null hypotheses using a convex decreasing density estimate on a vector of
p-values.
}
\usage{
convest(p, niter = 100, doplot = FALSE, doreport = FALSE)
}
\arguments{
  \item{p}{numeric vector of p-values, calculated using any test of your
    choice. Missing values are not allowed}
  \item{niter}{number of iterations to be used in fitting the convex,
    decreasing density for the p-values. Default is 100.}
  \item{doplot}{logical, should updated plots of fitted convex decreasing
    p-value density be produced at each iteration?
    Default is \code{FALSE}.}
  \item{doreport}{logical, should the estimated proportion be printed at each iteration?
  Default is \code{FALSE}.}
}
\details{
The proportion of true null hypotheses is often denoted pi0.
}
\value{
Numeric value in the interval [0,1] representing he estimated proportion of true null hypotheses, with \code{class} being \code{convest} and the \code{lfdr} attribute containing estimated local false discovery rates. 
}
\references{
Langaas, M., Ferkingstad, E., and Lindqvist, B. (2005).
Estimating the proportion of true null hypotheses, with application
to DNA microarray data.
\emph{Journal of the Royal Statistical Society Series} B, 67, 555-572.
Preprint at \url{http://www.math.ntnu.no/~mettela/}
}
\author{Long Qu slightly modified and \code{\link[limma:convest]{convest}} function by Egil Ferkingstad and Mette Langaas in \code{limma} package. }
\seealso{
See \link[limma:08.Tests]{08.Tests} for other functions for producing or interpretting p-values.
}

\examples{
set.seed(9992722)
convest(runif(1e4)^1.5, niter=100,doplot=interactive())
}
\keyword{htest}

