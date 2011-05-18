\name{fitted.discTMix}
\Rdversion{1.1}
\alias{fitted.discTMix}
\alias{fitted.values.discTMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density evaluated at observed statistics for discTMix class
}
\description{Density evaluated at observed statistics for discTMix class
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{fitted}{discTMix}(object, ...)
\method{fitted.values}{discTMix}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{discTMix}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{Not used currently. 
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{numeric vector of densities at each observed statistic of \code{object}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Pawitan Y, Krishna Murthy KR, Michiels S, Ploner A (2005) Bias in the estimation of false discovery rate in microarray studies, Bioinformatics.
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
%\note{
%Functions for \code{nparncpp} are not yet implemented. 
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sparncpt}}, \code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{nparncpp}}, \code{\link{discTMix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% __ONLY ONE__ keyword per line
