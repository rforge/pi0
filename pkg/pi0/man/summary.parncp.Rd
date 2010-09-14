\name{summary.parncpt}
\Rdversion{1.1}
\alias{summary.parncpt}
\alias{print.parncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Print summary for an object of class parncp
}
\description{
        Print summary for an object of class \code{parncpt}. This includes the estimated proportion of true null hypotheses (\code{pi0}),
estimated mean of noncentrality parameters (\code{sd.ncp}), estimated standard deviation of noncentrality parameters (\code{sd.ncp}).
}
\usage{
\method{summary}{parncpt}(object,...)
\method{print}{parncpt}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,x}{ an object of class \code{parncpt}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
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

\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{sparncpt}}
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
