\name{dtn.mix}
\Rdversion{1.1}
\alias{dtn.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of noncental t-normal mixture
}
\description{Density of noncentral t-distribution, with noncentrality parameter (NCP) being normally distributed. This is a scaled noncentral t-density. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dtn.mix(t, df, mu.ncp, sd.ncp, log = FALSE, approximation = c("int2", 
        "saddlepoint", "laplace", "none"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{A numeric vector of quantiles
%%     ~~Describe \code{t} here~~
}
  \item{df}{A numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{mu.ncp}{A numeric vector of normal mean of NCP
%%     ~~Describe \code{mu.ncp} here~~
}
  \item{sd.ncp}{A numeric vector of normal SD of NCP
%%     ~~Describe \code{sd.ncp} here~~
}
  \item{log}{logical; if \code{TRUE}, log density is returned. 
%%     ~~Describe \code{log} here~~
}
  \item{approximation}{character; Method of approximation. \code{int2} computes exact denstiy for \code{int}eger \code{df} and polynomially \code{int}erpolate to non-integer degrees of freedom. 
  \code{saddlepoint} computes the saddle point approximation of the noncentral t-density. 
  \code{laplace} computes the laplacian approximation of the noncentral t-density. 
  \code{none} uses the (sort of) true noncentral t-density \code{\link{dt}} function. 
%%     ~~Describe \code{approximation} here~~
}
  \item{\dots}{other arguments passed to \code{\link{dt.int2}} or \code{\link{dt.sad}}. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{Mathematically, this is equivalent to \code{dt(t/s, df, mu.ncp/s)/s} where \code{s=sqrt(1+sd.ncp*sd.ncp)}. But the various approximations are usually sufficient for large problems where speed is more important than precision. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{numeric vector of densities
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Broda, Simon and Paolella, Marc S. (2007) Saddlepoint approximations for the doubly noncentral t distribution, Computational Statistics & Data Analysis, 51,6, 2907-2918.

Young, G.A. and Smith R.L. (2005) Essentials of statistical inference. Cambridge University Press. Cambridge, UK. 

Qu L, Nettleton D, Dekkers JCM. (2012) Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. Biometrics (in press).
}
\author{Long Qu 
%%  ~~who you are~~
}
\note{For normal-normal mixture, set \code{df=Inf}. When this is the case, \code{approximation} is ignored. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dt.sad}}, \code{\link{dt.int2}}, \code{\link{dt.lap}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
