\name{logLik.ncpest}
\Rdversion{1.1}
\alias{logLik.ncpest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{log likelihood from an object of class ncpest
}
\description{log likelihood from an object of class \code{ncpest}. This could be penalized liklihood in the case of a \code{nparncpt} object. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{logLik}{ncpest}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the object of class \code{ncpest}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{Extract the \code{logLik} component. This is used by \code{\link{AIC}}. The \code{df} is the estimated effetive number of parameters. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{an object of class \code{\link{logLik}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{Long Qu 
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{logLik}}, \code{\link{AIC}}, 
\code{\link{nparncpt}}, \code{\link{parncpt}}, \code{\link{sparncpt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{  }% __ONLY ONE__ keyword per line
