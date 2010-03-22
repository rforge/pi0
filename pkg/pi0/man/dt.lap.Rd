\name{dt.lap}
\Rdversion{1.1}
\alias{dt.lap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{laplacian approximation to noncentral t-distribution
%%  ~~function to do ... ~~
}
\description{laplacian approximation to noncentral t-distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.lap(x, df, ncp = 0, log = FALSE, normalize = c("central", "integral", "none"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector of quantiles
%%     ~~Describe \code{x} here~~
}
  \item{df}{numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{numeric vector of noncentrality parameters
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{logical; if \code{TRUE}, log densities are returned
%%     ~~Describe \code{log} here~~
}
  \item{normalize}{character. If \code{central}, the normalization is such that the approxiamtion is exact when \code{ncp=0}. 
  If \code{integral}, numerical integration is performed such that the density integrate to 1 (not implemented yet). 
  If \code{none}, no normalization is performed. 
%%     ~~Describe \code{normalize} here~~
}
  \item{\dots}{currently not used. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{
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
Young, G.A. and Smith R.L. 2005. Essentials of statistical inference. Cambridge University Press. Cambridge, UK. 
%% ~put references to the literature/web site here ~
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dt.int2}}, \code{\link{dt.sad}}, \code{\link{dt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(x, df, ncp=0, log=FALSE, normalize=c('central', 'integral','none'),...)
%{
%    x0=(sqrt(4*df*df+4*x*x*df+x*x*ncp*ncp)+x*ncp)/2/sqrt(x*x+df)
%    logC=df/2*log(df/2)-df*ncp*ncp/2/(x*x+df)-.5*log(pi/2)-lgamma(df/2)-(df+1)/2*log(x*x+df)
%    logNCT=logC+df*log(x0)-.5*(x0-ncp*x/sqrt(x*x+df))^2+.5*log(2*pi*x0*x0/(x0*x0+df))+
%            pnorm(0,x0,sqrt(1/(1+df/x0/x0)),lower=FALSE,log=TRUE)
%
%    normalization=match.arg(normalize)
%    if (normalization=='integral'){
%          warning('integral normalization not implemented yet; set normalize="central"')
%         normalization='central'
%    }
%    logNormalizer=if(normalization=='central') lgamma(df/2+1/2)+df/2-.5*log(2*pi)-df/2*log(df/2)-
%                                               pnorm(0,sqrt(df),1/sqrt(2),lower=FALSE,log=TRUE)
%                  else if (normalization=='none') 0 
%
%    logAns=logNormalizer+logNCT
%    if(log)return(logAns) else return(exp(logAns))
%
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar}
\keyword{ distribution }% __ONLY ONE__ keyword per line
