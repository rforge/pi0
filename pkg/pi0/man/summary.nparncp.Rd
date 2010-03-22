\name{summary.nparncp}
\Rdversion{1.1}
\alias{summary.nparncp}
\alias{print.nparncp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Print summary for an object of class nparncp
}
\description{
        Print summary for an object of class \code{nparncp}. This includes the estimated proportion of true null hypotheses (\code{pi0}),
estimated mean of noncentrality parameters (\code{sd.ncp}), estimated standard deviation of noncentrality parameters (\code{sd.ncp}), 
the estimated effective number of parameters (\code{enp}), and the tuning parameter (\code{lambda}) that controls the amount of smoothing. 
}
\usage{
\method{summary}{nparncp}(object, ...)
\method{print}{nparncp}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,x}{ an object of class \code{nparncp}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ the \code{\link{invisible}} \code{object} itself
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncp}}, \code{\link{nparncp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(object,...)
%{
%    print.nparncp(object,...)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ univar }% __ONLY ONE__ keyword per line
