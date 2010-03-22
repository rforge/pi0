\name{plot.parncp}
\Rdversion{1.1}
\alias{plot.parncp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot an object of class parncp, i.e., parametric estiamte of noncentrality parameters
%%  ~~function to do ... ~~
}
\description{Plot the histogram of observed t-statistics together with its fitted density estimate; also plotted is the estimated density of noncentrality parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{parncp}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{parncp}
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{ Left panel shows the density estimate of observed t-statistics, overlap with a histogram; right panel shows the estimated density of noncentrality parameters. Solid line is the actual mean of the estimate; dashed line is located at zero. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{the \code{\link{invisible}} \code{x} itself
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

\seealso{\code{\link{parncp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(x)
%{
%    x11(width=8, height=4)
%    op=par(mfrow=c(1,2))
%    hist(obj$data$tstat, pr=TRUE, br=min(c(max(c(20, length(obj$data$tstat)/100)), 200)), xlab='t',main='t-statistics')
%    ord=order(obj$data$tstat)
%    lines(obj$data$tstat[ord], fitted.parncp(x)[ord], col='red', lwd=3)
%    curve(dnorm(x, obj$mu.ncp, obj$sd.ncp), min(obj$data$tstat), max(obj$data$tstat), 500, xlab='delta', ylab='density',main='noncentrality parameters')
%    abline(v=c(0, obj$mu.ncp), lty=1:2)
%    par(op)
%    invisible(x)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ univar }% __ONLY ONE__ keyword per line
