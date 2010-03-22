\name{mTruncNorm}
\Rdversion{1.1}
\alias{mTruncNorm}
\alias{mTruncNorm.int2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moments of truncated normal distribution and the integral in the noncentral t-distribution
%%  ~~function to do ... ~~
}
\description{Compute the moments of truncated normal distribution and the integral that appears in the noncentral t-distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mTruncNorm(r = 1, mu = 0, sd = 1, lower = -Inf, upper = Inf, 
        approximation = c("int2", "laplace", "numerical"), integral.only = FALSE, ...)
mTruncNorm.int2(r = as.integer(1), mu = 0, sd = 1, lower = -Inf, upper = Inf, takeLog = TRUE, ndiv = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{the order of moments to be computed. It could be noninteger, but has to be nonnegative. This is also the degrees of freedom for the noncentral t-distribution.
%%     ~~Describe \code{r} here~~
}
  \item{mu}{mean of the normal distribution, before truncating. 
%%     ~~Describe \code{mu} here~~
}
  \item{sd}{SD of the normal distribution, before truncating. 
%%     ~~Describe \code{sd} here~~
}
  \item{lower}{lower truncation point
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{upper truncation point
%%     ~~Describe \code{upper} here~~
}
  \item{approximation}{Method of approximation. \code{int2} is exact for \emph{int}eger \code{r} and \emph{int}erpolate to noninteger \code{r}. 
\code{laplace} uses laplacian approximation. \code{numerical} uses nuemerical integration. 
%%     ~~Describe \code{approximation} here~~
}
\item{integral.only}{logical. If \code{TRUE}, only the integral in noncentral t-distribution is returned. Otherwise, it is normalized to be the rth moments of truncated normal distribution. 
%%     ~~Describe \code{integral.only} here~~
}
  \item{takeLog}{logical. If \code{TRUE} and \code{r} is not an integer, the polyomial interpolation will be on the log scale. But final result is on the original scale.
}
  \item{ndiv}{number of points with closes integer \code{r} to be used in polynomial interpolation. 
}
  \item{\dots}{other arguments passed to \code{mTruncNorm.int2}
%%     ~~Describe \code{\dots} here~~
}
}
\details{\code{mTruncNorm.int2} uses iterative relation over \code{r} to compute the integral iteratively starting from \code{r=0} and \code{r=1} whose analytic results are available. 
If \code{r} is not an integer, the nearest \code{ndiv} nonnegative integer \code{r} will be used to do divided difference polynomial interpolation. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{numeric vector. If \code{integral.only} is \code{TRUE}, this is the integral in the noncentral t-density; otherwise this is the rth moments of truncated normal distribution. 
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

\seealso{\code{\link{dt}}, \code{\link{pt}}, \code{\link{dt.int2}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(r=1, mu=0, sd=1, lower=-Inf, upper=Inf, approximation=c('int2','laplace', 'numerical'),integral.only=FALSE,...)
%{   ## return E(X^r|lower<X<upper), where X~Normal(mu, sd), integral.only=FALSE (default); 
%    ## return int_lower^upper X^r exp(-0.5*(X-mu)^2/sd^2) dx, if integral.only=TRUE
%    if(any(lower>upper)) stop("boundary inconsistent")
%    approximation=match.arg(approximation)
%    if(integral.only) fact=1 else fact=1 /  (pnorm(upper,mu,sd)-pnorm(lower,mu,sd))/sqrt(2*pi)/sd 
%
%    fact *
%    if(approximation=='laplace'){
%        x0=(sqrt(mu*mu+4*sd*sd*r)+mu)/2
%        x0^r*exp(-(x0-mu)*(x0-mu)/2/sd/sd)*
%            sqrt(2*pi/(r/x0/x0+1/sd/sd))*
%            (pnorm(upper,x0,1/sqrt(r/x0+1/sd/sd))- pnorm(lower,x0,1/sqrt(r/x0+1/sd/sd))) 
%    }else if (approximation=='int2') {
%        mTruncNorm.int2(r,mu,sd,lower,upper,...)
%    }else if (approximation=='numerical') {
%        obj=function(x,r,mu,sd)x^r*exp(-(x-mu)*(x-mu)/2/sd/sd) ## not working for negative values
%        n=max(c(length(r), length(mu), length(sd), length(lower), length(upper)))
%         mu=rep(as.double(mu), length=n); sd=rep(as.double(sd),length=n); 
%        lower=rep(as.double(lower), length=n);     upper=rep(as.double(upper),length=n)
%       ans=numeric(n)
%        for(i in 1:n) ans[i]=integrate(obj, lower[i], upper[1],r=r[i],mu=mu[i],sd=sd[i])$value  
%        ans
%    }
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
