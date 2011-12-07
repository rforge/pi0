\name{pdf.dist}
\alias{pdf.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Distance between densities }
\description{
        Compute the distance between two density functions
}
\usage{
pdf.dist(f1, f2, method = c("Hellinger", "abdif"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f1,f2}{ Two functions of one argument, both of which are densities defined over the whole real line. }
  \item{method}{ character; specifying the definition of distance. Current choices are \code{Hellinger} and \code{abdif}, 
  the latter of which is the integrated absolute differences between the two function. }
}
\details{
        Numerical integration is performed from \code{-Inf} to \code{Inf}. Hence, the two functions must be able to accept 
        arguments over the whole real line. 
}
\value{
        a numeric scalar which is the computed distance, or \code{NA_real_} if any problem occurs. 
}
%\references{
%}
\author{Long Qu }
%\note{
%
%
%}
%\seealso{  }
\examples{
# Hellinger distance betwee standard normal and log-normal
pdf.dist(dnorm, dlnorm, 'Hell') # 0.5981035

# absolute difference between standard normal and standard cauchy
f2=function(x)dt(x,1)
pdf.dist(dnorm, f2, 'abd')  #[1] 0.5023312
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
