\name{bw.mse.pdf.asym}
\Rdversion{1.1}
\alias{bw.mse.pdf.asym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find asymptotically optimal bandwidth for KDE by minimizing MSE at the first data point
}
\description{
Find asymptotically optimal bandwidth for Gaussian kernel density estimation by minimizing mean squared error at the first data point. 
}
\usage{
bw.mse.pdf.asym(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector whose density is to be estimated. 
}
  \item{iter.max}{The max number of iterations. For an iterative procedure, set this to be  greater than 1L.}
  \item{eps}{A small positive number such that iterations are terminated when the difference between successive bandwidths is below \code{eps}.
}
  \item{start.bw}{A function, a character, or a positive number. If being a function or a character, the corresponding function will be called with \code{x} as the only argument to get the initial bandwidth. When it is a number, it will be used as the starting bandwidth directly.
}
  \item{verbose}{Logical, when \code{TRUE}, the bandwidth for each iteration will be printed. 
}
}
\details{
The procedure starts from the starting bandwidth specified by \code{start.bw}. The new bandwidth will be computed as 
\deqn{
\left[\frac{Rf(x_1)}{N\lbrace f''(x_1) \rbrace^2} \right]^{1/5}
}{
( {Rf(x_1)} / [N{f''(x_1)}^2] )^{1/5}
}
where \code{R} is approximately 0.2821, \code{N} is the sample size and \code{f} is the density estimate using the current bandwidth. 
}
\value{
A numeric positive scalar, giving the final bandwidth. 
}
\references{
Parzen, E. (1962). On Estimation of a Probability Density Function and Mode. The Annals of Mathematical Statistics 33: 1065--1076.

Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization. Wiley, New York.

Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall, London. 
}
\author{Long Qu
}
\seealso{
\code{\link{density}}, \code{\link[ks]{kdde}}
}
\examples{
set.seed(23490)
x=rnorm(100)
bw.mse.pdf.asym(x)
\dontrun{
plot(density(x, bw.mse.pdf.asym(x)), ylim=c(0,dnorm(0)))
abline(v=x[1L])
curve(dnorm, lty=3, col='grey', add=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
