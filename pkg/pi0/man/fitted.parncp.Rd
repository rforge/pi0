\name{fitted.parncp}
\Rdversion{1.1}
\alias{fitted.parncp}
\alias{fitted.values.parncp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density evaluated at observed t-statistics for parncp class
%%  ~~function to do ... ~~
}
\description{Density evaluated at observed t-statistics for \code{parncp} class
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{fitted}{parncp}(object, ...)
\method{fitted.values}{parncp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{parncp}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{other arguments passed to \code{dtn.mix}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{numeric vector of densities at each t-statistic of \code{object}
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

\seealso{\code{\link{fitted.nparncp}}, \code{\link{parncp}}, \code{\link{nparncp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(object, ...)
%{
%    obj$pi0*dt(obj$data$tstat, obj$data$df)+(1-obj$pi0)*dtn.mix(obj$data$tstat, obj$data$df,obj$mu.ncp,obj$sd.ncp,FALSE,...)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
