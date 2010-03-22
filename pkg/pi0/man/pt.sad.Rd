\name{Saddle Point Approximation of noncentral t-distribution}
\Rdversion{1.1}
\alias{pt.sad}
\alias{dt.sad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Saddle Point Approximation of noncentral t-distribution
%%  ~~function to do ... ~~
}
\description{ Density and cumulative distribution function of noncentral t-distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.sad(x, df, ncp = 0, log = FALSE, normalize = c("approximate", "derivative", "integral", "none"), epsilon = 1e-04)
pt.sad(q, df, ncp = 0, log = FALSE, epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{ numeric vector of quantiles
%%     ~~Describe \code{q} here~~
}
  \item{df}{ numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{ numeric vector of noncentrality parameter
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{ logical; whether \code{log} should be taken. 
%%     ~~Describe \code{log} here~~
}
  \item{epsilon}{ a small numeric scalar; if the difference between \code{q} and \code{ncp} is closer than this, results will be computed differently. 
%%     ~~Describe \code{epsilon} here~~
}
  \item{normalize}{ the way to normalize the approximate density so that it is closer to a true density. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ \code{dt.sad} returns density; \code{pt.sad} returns the probability. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Broda, Simon and Paolella, Marc S. 2007. Saddlepoint approximations for the doubly noncentral t distribution, Computational Statistics & Data Analysis, 51,6, 2907-2918.
%% ~put references to the literature/web site here ~
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{dt.int2}}, \code{\link{dt.lap}}, \code{\link{dt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(q, df, ncp=0, log=FALSE,epsilon=1e-4){
%    if(missing(df)) stop('Error: the number of degrees of freedom must be supplied!')
%
%    N=max(c(length(q),length(df),length(ncp)))
%    y1=rep(q,length=N)
%    mu=ncp
%    n=n
%
%    y2.hat=(mu*y1+sqrt(4*n*(y1*y1+n)+mu*mu*y1*y1))/2/(y1*y1+n)
%    t1.hat=-mu+y1*y2.hat
%    t2.hat=-y1*t1.hat/2/n/y2.hat
%    w=sqrt(-mu*t1.hat-2*n*log(y2.hat))*ifelse(y1>mu,1,-1)
%    u=sqrt(mu*y1*y2.hat/2/n+1)/y2.hat
%    d=1/t1.hat/y2.hat
%    alpha=mu
%
%    ans=pnorm(w)+dnorm(w)*(1/w-d/u)
%
%    idx=abs(y1-alpha)<epsilon
%    if(any(idx)){
%        lim.F=1/2-1/6/sqrt(n*pi*(mu*mu+2*n)^3)*mu*(2*mu*mu+3*n)
%        n2=n*n
%        n3=n2*n
%        n4=n2*n2
%        n5=n2*n3
%        o1=6*n*n-n
%        o2=144*n*n-33*n
%        o3=4*n*n-n
%        mu2=mu*mu
%        mu4=mu2*mu2
%        mu6=mu2*mu4
%        lim.F.prime=(2*mu6*n*o1+12*mu4*o1*n2+mu2*n3*o2+24*n4*o3)/(12*sqrt(pi*n5*(mu2+2*n)^7))
%
%        ans[idx]=(lim.F+(y1-alpha)*lim.F.prime)[idx]
%    }
%    if(log)return(log(ans)) else return(ans)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
%\keyword{ univar }% __ONLY ONE__ keyword per line
