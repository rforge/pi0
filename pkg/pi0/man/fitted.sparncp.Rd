\name{fitted.sparncpt}
\Rdversion{1.1}
\alias{fitted.sparncpt}
\alias{fitted.values.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density evaluated at observed t-statistics for sparncpt class
}
\description{Density evaluated at observed t-statistics for \code{sparncpt} class
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{fitted}{sparncpt}(object, ...)
\method{fitted.values}{sparncpt}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{sparncpt}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{other arguments passed to \code{dtn.mix}
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{numeric vector of densities at each t-statistic of \code{object}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)
}
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
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
