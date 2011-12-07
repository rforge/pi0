\name{lastbin}
\alias{lastbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Histogram estimator of p-value density evaluated at 1 }
\description{
This function reports the density estimate of the right most bin of histogram of p-values.
}
\usage{
lastbin(p, bw = 0.2, trunc = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ a numeric vector the p-values }
  \item{bw}{ numeric, the bin width of histogram of p-values }
  \item{trunc}{ logical, indicating if the resulting estimate should be truncated to within [0,1]. }
}
\details{
This is a very fast and cheap estimate of p-value density at one, with a slight positive bias, 
because it is only an unbiased estimate of density at the bin center, not the right bin edge. 
But this is usually ignorable. The function is defined as:\cr
\code{function(p,bw=.2,trunc=TRUE)
    if(trunc)max(min(1,mean(p>=1-bw)/bw),0) else mean(p>=1-bw)/bw
}.
}
\value{
	a single numeric value as the estimate. 
}

\author{ Long Qu  }
\seealso{\code{\link[qvalue:qvalue]{qvalue}}, \code{\link{histf1}}}
\examples{
set.seed(9992722)
lastbin(runif(5e5)^1.5) ##  [1] 0.69511
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
