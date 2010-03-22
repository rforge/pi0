\name{dt.int2}
\Rdversion{1.1}
\alias{dt.int2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{INTerpolation of INTeger degrees of freedom noncentral t-density
%%  ~~function to do ... ~~
}
\description{This funciton evaluate the noncentral t-density using an iterative procedure for integer degrees of freedom. This is much faster than two calls to the \code{pt} approach. 
For non-integer degrees of freedom, the polynomial interpolation is used to approximate the density. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.int2(x, df, ncp, log = FALSE, ndiv = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector of quantiles
%%     ~~Describe \code{x} here~~
}
  \item{df}{numeric vector degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{numeric vector of noncentrality parameters
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{logical; if \code{TRUE}, log densities are returned. 
%%     ~~Describe \code{log} here~~
}
  \item{ndiv}{numeric; the number of points used for polynomial interpolation
%%     ~~Describe \code{ndiv} here~~
}
}
\details{This function uses the iterative relation for the integral in the noncentral t-density. It starts with \code{df=0} and \code{df=1}, 
and then iteratively computes the integral for larger \code{df}. For non-integer \code{df}, it uses \code{ndiv} nearest points to perform a 
divided difference polynomial interpolation. For integer \code{df}, this function is about 2 to 3 times faster than \code{dt} function. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{numeric vector of densities. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dtn.mix}}, \code{\link{dt}}, \code{\link{dt.lap}}, \code{\link{dt.sad}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(x, df, ncp, log=FALSE, ndiv=8 ) ## pretty fast computation of noncentral t density when df is integer
%{   ## when df is integer, this is exact for noncentral t density;
%    ## when df is fractional, this is divided difference polynomial interpolation using ndiv points with nearest integer dfs
%    if (missing(ncp)) 
%        return(.Internal(dt(x, df, log)))
%
%    n=max(c(length(x),length(df),length(ncp)))
%    x=rep(as.double(x),length=n); df=rep(df,length=n); ncp=rep(as.double(ncp),length=n)
%
%    integer.df.idx=(df==round(df))
%    n.int=as.integer(sum(integer.df.idx))
%
%    df.half=df/2
%    tsq.df=x*x+df
%
%    logC=df.half*log(df.half)-0.5*log(pi/2)-lgamma(df.half)-(df.half+.5)*log(tsq.df)-df.half/tsq.df*ncp*ncp
%    mus=as.double(x/sqrt(tsq.df)*ncp)
%
%    ints=numeric(n); ints+0
%
%    if(n.int>0)   ints[integer.df.idx]=.C('intTruncNormVec', n.int, 
%                                               as.integer(df[integer.df.idx]), 
%                                               mus[integer.df.idx],  rep(as.double(1.0), n.int), 
%                                               numeric(n.int), rep(Inf,n.int), ans=ints[integer.df.idx], NAOK=TRUE)$ans
%    if(n-n.int>0) ints[!integer.df.idx]=.C('fracTruncNormVec', n=n-n.int, r=as.double(df)[!integer.df.idx], 
%                            mu=mus[!integer.df.idx], sd=rep(as.double(1.0),n-n.int), 
%                            low=rep(as.double(0.0), n-n.int), upp=rep(Inf, n-n.int),
%                            ans=ints[!integer.df.idx], ndiv=as.integer(ndiv), takeLog=as.integer(1), NAOK=TRUE)$ans
%
%    if(log) { logC+log(pmax(0,ints))
%    }else exp(logC)*pmax(0,ints)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
