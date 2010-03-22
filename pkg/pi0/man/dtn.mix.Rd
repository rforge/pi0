\name{dtn.mix}
\Rdversion{1.1}
\alias{dtn.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of noncental t-normal mixture
%%  ~~function to do ... ~~
}
\description{Density of noncentral t-distribution, with nencentrality parameter (NCP) being normally distributed. This is a scaled noncentral t-density. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dtn.mix(t, df, mu.ncp, sd.ncp, log = FALSE, approximation = c("int2", "saddlepoint", "laplace", "none"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{numeric vector of quantiles
%%     ~~Describe \code{t} here~~
}
  \item{df}{numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{mu.ncp}{numeric vector of normal mean of NCP
%%     ~~Describe \code{mu.ncp} here~~
}
  \item{sd.ncp}{numeric vector of normal SD of NCP
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
\details{Mathematically, this is equivalent to \code{dt(t/s, df, mu.ncp/s)/s} where \code{s=sqrt(1+sd.ncp*sd.ncp)}. But various approximation is usually sufficient for large problems where speed is more important than precision. 
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
Broda, Simon and Paolella, Marc S. 2007. Saddlepoint approximations for the doubly noncentral t distribution, Computational Statistics & Data Analysis, 51,6, 2907-2918.
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

\seealso{\code{\link{dt.sad}}, \code{\link{dt.int2}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(t,df,mu.ncp, sd.ncp, log=FALSE, approximation=c('int2','saddlepoint','laplace','none'),...)
%{
%    approximation=match.arg(approximation)
%    if(approximation=='none'){
%        scale.fact=sqrt(1+sd.ncp*sd.ncp)
%        ncp=mu.ncp/scale.fact
%        return( if(log){
%                    dt(t/scale.fact, df, ncp, log=TRUE)-log(scale.fact)
%                }else{
%                    dt(t/scale.fact, df, ncp)/(scale.fact)
%                }
%        )
%    }else if (approximation=='int2'){
%        scale.fact=sqrt(1+sd.ncp*sd.ncp)
%        ncp=mu.ncp/scale.fact
%        return( if(log){
%                    dt.int2(t/scale.fact, df, ncp, log=TRUE,...)-log(scale.fact)
%                }else{
%                    dt.int2(t/scale.fact, df, ncp,...)/(scale.fact)
%                }
%        )
%    }else if (approximation=='saddlepoint'){
%        scale.fact=sqrt(1+sd.ncp*sd.ncp)
%        ncp=mu.ncp/scale.fact
%        return( if(log){
%                    dt.sad(t/scale.fact, df, ncp, log=TRUE,...)-log(scale.fact)
%                }else{
%                    dt.sad(t/scale.fact, df, ncp,...)/(scale.fact)
%                }
%        )
%    }else if (approximation=='laplace'){
%        denom=(1+sd.ncp*sd.ncp)*df+t*t
%        u0=mu.ncp*t*sqrt(df+t*t)/denom
%        g0sq=(df+t*t)*(1+sd.ncp*sd.ncp)/denom
%        g0=sqrt(g0sq)
%        x0=(sqrt(4*g0sq*df+u0*u0)+u0)/2
%
%        exponent=-mu.ncp*mu.ncp*df/2/denom + df*log(x0) -(x0-u0)*(x0-u0)/2/g0sq
%
%        norm.prob=pnorm(0, x0, x0*g0/sqrt(g0sq*df+x0*x0), lower=FALSE, log=TRUE)
%        norm.prob0=pnorm(0,g0*sqrt(df),g0/sqrt(2),lower=FALSE,log=TRUE)
%
%
%        ans=(exponent
%            +log(x0)+df/2-df/2*log(df+t*t)
%            -lbeta(df/2,.5)
%            -.5*log(.5)-.5*log((1+sd.ncp*sd.ncp)*df+t*t)-.5*log(g0sq*df+x0*x0)
%            +norm.prob - norm.prob0
%        )
%        return (    if(log) ans else exp(ans)   )
%    }
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
