\name{plot.nparncpt}
\Rdversion{1.1}
\alias{plot.nparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot an object of class nparncpt, i.e., nonparametric estiamte of noncentrality parameters
%%  ~~function to do ... ~~
}
\description{Plot the Network information criterion (NIC), effective number of parameters (ENP), and estimated proportion (pi0) of true null hypothesis
for different choice of tuning parameters; also plot the estimated density of noncentrality parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{nparncpt}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{nparncpt}
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{currently not used. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{For NIC, only values within one s.e. of the minimum are shown. The solid line on NIC, ENP and pi0 shows the final tuning parameter, i.e., the one that minimizes NIC. 
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
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{nparncpt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(x,...)
%{
%    x11(width=7,heigh=7)
%    op=par(mfrow=c(2,2))
%    attach(x)
%    n.lambda=length(all.lambdas)
%    i.1se=tail(which(all.nics<=all.nics[i.final]+all.nic.sd[i.final] & 1:n.lambda>=i.final),1) ## one se rule
%    plot(log10(all.lambdas), (all.nics), 
%        ylim=range(c(all.nics+all.nic.sd, all.nics-all.nic.sd)[c(1:min(i.1se,i.final+1,n.lambda),n.lambda+1:min(i.1se,i.final+1,n.lambda))]), 
%        xlab='log10(lambda)', ylab='NIC')
%    for(i in 1:n.lambda) lines(rep(log10(all.lambdas)[i],2), all.nics[i]+c(-1,1)*all.nic.sd[i], lwd=3)
%    abline(v=log10(all.lambdas[i.final]), col=2);
%    abline(h=all.nics[i.final]+c(-1,1)*all.nic.sd[i.final], col=4)
%#    abline(v=log10(all.lambdas[i.1se]), col=2,lty=2)
%
%#    plot(log10(all.lambdas), all.enps); abline(v=log10(all.lambdas[c(i.final,i.1se)]), xlab='log10(lambda)', ylab='ENP', lty=1:2)
%#    plot(log10(all.lambdas), all.pi0s); abline(v=log10(all.lambdas[c(i.final,i.1se)]), xlab='log10(lambda)', ylab='pi0', lty=1:2)
%    plot(log10(all.lambdas), all.enps,xlab='log10(lambda)',ylab='effective # parameters'); 
%        abline(v=log10(all.lambdas[c(i.final)]), xlab='log10(lambda)', ylab='ENP', lty=1)
%    plot(log10(all.lambdas), all.pi0s, xlab='log10(lambda)', ylab='pi0'); 
%        abline(v=log10(all.lambdas[c(i.final)]), xlab='log10(lambda)', ylab='pi0', lty=1)
%
%    d.ncp=function(x, betas)        # p in the paper
%    {   ## depends on mus, sigs
%        x=outer(x, all.mus, '-')
%        x=sweep(x, 2, all.sigs, '/')
%        d=sweep(dnorm(x),2,all.sigs,'/')
%        drop(d\%*\%betas)
%    }
%    curve(d.ncp(x, beta), min(data$tstat),max(data$tstat),100,col=4,lwd=2, xlab='delta',ylab='density')
%    detach(x)
%    par(op)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot}
\keyword{ univar }% __ONLY ONE__ keyword per line
