\name{summary.sparncp}
\alias{summary.sparncp}
\alias{print.sparncp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{        Print summary for an object of class sparncp
 }
\description{
        Print summary for an object of class \code{sparncp}. This includes the estimated proportion of true null hypotheses (\code{pi0}),
estimated mean of noncentrality parameters (\code{sd.ncp}), estimated standard deviation of noncentrality parameters (\code{sd.ncp}), and 
effective number of parameters (\code{enp}).
}
\usage{
\method{summary}{sparncp}(object, ...)
\method{print}{sparncp}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,x}{ an object of class \code{sparncp}
  \item{\dots}{ currently not used }
}
\details{
}
\value{ the \code{\link{invisible}} \code{object} itself
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncp}}, \code{\link{nparncp}}, \code{\link{sparncp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(object)
%{
%    cat("pi0 (proportion of null hypotheses) =", obj$pi0, fill=TRUE)
%    cat("mu.ncp (mean of noncentrality parameters) =", obj$mu.ncp, fill=TRUE)
%    cat("sd.ncp (SD of noncentrality parameters) =", obj$sd.ncp, fill=TRUE)
%    invisible(object)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ univar }% __ONLY ONE__ keyword per line
