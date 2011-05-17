\name{CBUM}
\alias{CBUM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ (Censored) Beta-Uniform mixture for p-values }
\description{
This function implements the method of Markitsis and Lai (2010).  
}
\usage{
CBUM(p, start.pi0=0.5, thresh.censor=0.05, eps=1e-5, niter=Inf, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ a numeric vector the p-values }
  \item{start.pi0}{ numeric scalar, starting value of pi0 for EM algorithm. }
  \item{thresh.censor}{numeric scalar, the threshold of censoring. If \code{isTRUE(lambda<min(p))}, this is equivalent to the BUM method of Pounds and MOrris (2003). }
  \item{eps}{numeric scalar, maximum tolerable absolute difference of parameter estimates for successive iterations in the EM algorithm. }
  \item{niter}{numeric scalar, maximum number of EM iterations.}
  \item{verbose}{logical scalar, indicating whether excessive outputs will be printed during EM algorithm.}
}
\details{
This function is an improved version of the \code{CBpi0} function available at http://home.gwu.edu/~ylai/research/CBpi0/CBpi0.txt, which
implements the censored (1-paramater beta)-uniform mixture model to a large number of p-values. 
}
\value{
	A numeric scalar, being the proportion pi0 of true null hypotheses. 
        The result has a class 'CBUM', with the following attributes:
        \item{converged}{logical, convergence status.}
        \item{iter}{numeric, number of iterations.}
        \item{call}{the \code{match.call()} result.}
        \item{alpha}{estimated alpha parameter for the beta component.}
        \item{lfdr}{numeric vector estimated local false discovery rates, if \code{thresh.censor<min(p)}; NULL, otherwise}
        \item{thresh.censor}{the censoring threshold for p-values.}        
}
\references{
Anastasios Markitsis and Yinglei Lai (2010) A censored beta mixture model for the estimation of the proportion of non-differentially expressed genes. Bioinformatics 26(5):640-646.
}
\author{ Long Qu \email{longor@iastate.edu} modified the code from Markitsis and Lai (2010). }
\seealso{\code{\link[qvalue:qvalue]{qvalue}}, \code{\link{histf1}}}
\examples{
set.seed(9992722)
p=c(runif(7000), rbeta(3000, .8, 2.2))
CBUM(p, eps=1e-5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
