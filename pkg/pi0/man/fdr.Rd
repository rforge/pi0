\name{fdr}
\alias{fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ False Discovery Rate (FDR) Estimation Based on a Given Estimate of Pi0 }
\description{
	This function estimate the \code{\link[qvalue]{qvalue}} based on p-values and an estimate of the 
	proportion of true null hypotheses, pi0. 
}
\usage{
fdr(p, pi0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ a numeric vector of p-values }
  \item{pi0}{ numeric, a given estimate of the proportion of true null hypotheses, \eqn{\pi_0}{pi0}, truncated to [0,1]. 
		The default \code{pi0=1} is the conservative \cite{Benjamini and Hochberg (1995)} version. }
}
\details{
	The estimation of q-value/FDR is the simple and quick plug-in method: 
	\deqn{q_i =\min_{i \le j \le G} G*\pi_0 * p_{(j)}/j }{q_i = min_{i<=j<=G} G*pi0*p_(j)/j}. 
}
\value{

a numeric vector of the same length as \code{p}, giving the estimated q-values corresponding to each p-value. 
}

\note{
This implementation avoids explicit loops and is much faster when the number of p-values are very large. 
}

\references{
Benjamini, Y., Hochberg, Y. (1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society, Series B (Methodological) 57, 289-300 \cr
Storey, J.D.,  Tibshirani, R. (2003) Statistical significance for genomewide studies. The Proceedings of the National Academy of Sciences 100, 9440-9445
}
\author{ Long Qu \email{lqu@biostatsolutions.com} }

\seealso{ \code{\link[qvalue:qvalue]{qvalue}}, \code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}}}
\examples{
set.seed(9992722)
pvals=runif(5e5)^1.5 ## simulate some fake 'p-values'
library(qvalue)
qvalObj=qvalue(pvals)  ## warning: this may be slow! 
fdrObj=fdr(pvals,qvalObj$pi0)
all.equal(fdrObj,qvalObj$qval) ## should be TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
