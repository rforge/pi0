\name{fitted}
\Rdversion{1.1}
\alias{fitted.nparncpt}
%\alias{fitted.values.nparncpt}
\alias{fitted.nparncpp}
%\alias{fitted.values.nparncpp}
\alias{fitted.parncpt}
\alias{fitted.parncpt2}
%\alias{fitted.values.parncpt}
\alias{fitted.sparncpt}
%\alias{fitted.values.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density evaluated at observed statistics for ncpest class
}
\description{Density evaluated at observed statistics for ncpest class
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{fitted}{nparncpt}(object, ...)
%\method{fitted.values}{nparncpt}(object, ...)
\method{fitted}{parncpt}(object, ...)
\method{fitted}{parncpt2}(object, ...)
%\method{fitted.values}{parncpt}(object, ...)
\method{fitted}{sparncpt}(object, ...)
%\method{fitted.values}{sparncpt}(object, ...)
\method{fitted}{nparncpp}(object, ...)
%\method{fitted.values}{nparncpp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{nparncpt}, \code{parncpt}, \code{sparncpt} or \code{nparncpp}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{other arguments passed to \code{dtn.mix}
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
Qu L, Nettleton D, Dekkers JCM. (2012) Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of \eqn{t}-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. Biometrics. 68. 1178-1187.
}
\author{Long Qu 
%%  ~~who you are~~
}
\note{
Functions for \code{nparncpp} are not yet implemented. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sparncpt}}, \code{\link{parncpt}}, \code{\link{parncpt2}}, \code{\link{nparncpt}}, \code{\link{nparncpp}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }% __ONLY ONE__ keyword per line
