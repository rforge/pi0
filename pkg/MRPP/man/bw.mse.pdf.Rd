\name{bw.mse.f}
\Rdversion{1.1}
\alias{bw.mse.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bw.mse.f(z, mc.n = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{mc.n}{
%%     ~~Describe \code{mc.n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(z, mc.n=500)
{
    K=mc.n
    y=qnorm(1:K/(K+1))
#    Phi.y=pnorm(y)
    phi.y=dnorm(y)
#    dphi.y=-y*phi.y
    t=z[1]
    h0=bw.nrd(z)
    B=length(z)
    t.z.h0=(t-z)/h0
    y.h0=y/h0
    sum.phi=sum(dnorm(t.z.h0))
    inner=sum.y.Phi.inner=sum.y.phi.Phi.inner=current.h=Inf    ## to be modified by obj or der.obj
    obj=function(h){
        if(current.h!=h){
            current.h<<-h
            inner<<-outer(t.z.h0, h*y.h0, '-')
                Phi.inner=pnorm(inner)
            sum.y.Phi.inner<<-sum(sweep(Phi.inner, 2, y, '*'))
            sum.y.phi.Phi.inner<<-sum(sweep(Phi.inner, 2, y*phi.y, '*'))
        }
      -2/B/B/h/h/K*sum.y.phi.Phi.inner+
      (B-1)/B/B/B/h/h/K/K*sum.y.Phi.inner*sum.y.Phi.inner+
      2/B/B/h/h0/K*sum.phi*sum.y.Phi.inner
    }
    der.obj=function(h){
        if(current.h!=h){
            current.h<<-h
            inner<<-outer(t.z.h0, h*y.h0, '-')
                Phi.inner=pnorm(inner)
            sum.y.Phi.inner<<-sum(sweep(Phi.inner, 2, y, '*'))            # not used in der.obj, but keep them up-to-date
            sum.y.phi.Phi.inner<<-sum(sweep(Phi.inner, 2, y*phi.y, '*'))         # not used in der.obj, but keep them up-to-date
        }
        phi.inner=dnorm(inner)
        sum.y.yh0.phi.inner=sum(sweep(phi.inner,2,y.h0*y,'*'))

        
     2/B/B/K*(sum(sweep(phi.inner,2,y.h0*y*phi.y,'*'))/(h*h)+2*sum.y.phi.Phi.inner/(h*h*h))-
     (B-1)/B/B/B/K/K*(2*sum.y.Phi.inner/h)*(sum.y.yh0.phi.inner/h+
                                            sum.y.Phi.inner/(h*h))-
     2*sum.phi/B/B/K/h0*(sum.y.yh0.phi.inner/h+
                         sum.y.Phi.inner/(h*h))
    }
    ofit=optim(bw.mse.f.asym(z,iter=FALSE), obj, der.obj, method='BFGS',
                control=list(fnscale=1e-6#,factr=1e3,parscale=1e-5
                ))
    if(ofit$par>0 && ofit$convergence==0){return(ofit$par)}
#    ofit=optim(bw.mse.f.asym(z,iter=FALSE), obj, der.obj, method='L-BFGS-B', lower=0, 
#                control=list(fnscale=1e-6#,factr=1e3,parscale=1e-5
#                ))
#    if(ofit$par>0 && ofit$convergence==0){return(ofit$par)}
    ofit=constrOptim(bw.mse.f.asym(z,iter=FALSE), obj, der.obj, matrix(1),0, 
                control=list(fnscale=1e-6))
    if(ofit$par>0 && ofit$convergence==0){return(ofit$par)}
    ans=ofit$par
    attr(ans,'optim')=ofit
    warning('nonconvergence in finding optimal bandwidth')
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
