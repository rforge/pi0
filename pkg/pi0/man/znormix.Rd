\name{znormix}
\alias{znormix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Normal-mixture based estimation of LFDR and pi0 }
\description{
This function implements the method of McLachLan, Bean and Jones (2006). 
}
\usage{
znormix(p, theoretical.null=TRUE, start.pi0, eps=1e-5, niter=Inf, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ a numeric vector the p-values }
  \item{theoretical.null}{ logical scalar, indicating whether theoretical N(0,1) null distribution is assumed for z-scores. }
  \item{start.pi0}{ optional numeric scalar, starting value of pi0 for EM algorithm; if missing, \code{\link[qvalue:qvalue]{qvalue}} will be called with default arguments to get this starting value. }
  \item{eps}{numeric scalar, maximum tolerable absolute difference of parameter estimates for successive iterations in the EM algorithm. }
  \item{niter}{numeric scalar, maximum number of EM iterations.}
  \item{verbose}{logical scalar, indicating whether excessive outputs will be printed during EM algorithm.}
}
\details{
A two-component normal mixture model is fit thru EM algorithm on the z-scores, where \code{z=qnorm(1-p)}. 
}
\note{
There are two small differences with McLachlan, Bean and Jones (2006): 
        \itemize{
        \item{If \code{start.pi0} is missing, it is estimated by the q-value smoother method implimented in \code{\link[qvalue:qvalue]{qvalue}}.}
        \item{For the empirical null case, a call to \code{quantile(z, start.pi0)} is used as the threshold to determine the initial component assignment, when choosing starting values. }
        }
}
\value{
	A length 5 numeric named vector of estimated parameters, with class 'znormix' and attributes
        
        \item{theoretical.null}{the same as input.}
        \item{converged}{logical, convergence status.}
        \item{iter}{numeric, number of iterations.}
        \item{call}{the \code{match.call()} result.}
        \item{lfdr}{numeric vector of local false discovery rates, with order being the same as the input p-values.}
        \item{fdr}{numeric vector of false discovery rates, with order being the same as the input p-values.}
        
}
\references{
G.J. McLachlan, R.W. Bean and L. Ben-Tovim Jones. (2006) A Simple implementation of a normal mixture approach to differential gene expression in multiclass microarrays. Bioinformatics, 22(13):1608-1615. 
}
\author{ Long Qu }
\seealso{\code{\link[qvalue:qvalue]{qvalue}}, \code{\link{histf1}}}
\examples{
set.seed(9992722)
p=1-pnorm(c(rnorm(7000),rnorm(3000,1)))
znormix(p, eps=1e-5)['pi0'] ##  0.6759555
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
