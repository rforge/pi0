\name{coef.ncpest}
\Rdversion{1.1}
\alias{coef.ncpest}
\alias{coefficients.ncpest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Extract estimated parameters
}
\description{
        Extract estimated parameters from object of \code{ncpect} class
}
\usage{
\method{coef}{ncpest}(object,...)
\method{coefficients}{ncpest}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of \code{ncpest} class
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{ a numeric vector of the estimated parameters
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

\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{sparncpt}}, \code{\link{nparncpp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function(object)
%{
%    obj$par
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ estimation }% __ONLY ONE__ keyword per line
