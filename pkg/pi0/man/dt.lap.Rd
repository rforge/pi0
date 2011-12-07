\name{dt.lap}
\Rdversion{1.1}
\alias{dt.lap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Laplace approximation to the noncentral t-density
}
\description{Laplace approximation to noncentral t-density
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.lap(x, df, ncp = 0, log = FALSE, normalize = c("central", 
      "integral", "none"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector of quantiles
%%     ~~Describe \code{x} here~~
}
  \item{df}{A numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{A numeric vector of noncentrality parameters
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{logical; if \code{TRUE}, log densities are returned
%%     ~~Describe \code{log} here~~
}
  \item{normalize}{character. If \code{central}, the normalization is such that the approxiamtion is exact when \code{ncp=0}. 
  If \code{integral}, numerical integration is performed such that the density integrate to 1 (not implemented yet). 
  If \code{none}, no normalization is performed. 
%%     ~~Describe \code{normalize} here~~
}
  \item{\dots}{currently not used. 
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{numeric vector of densities
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Young, G.A. and Smith R.L. 2005. Essentials of statistical inference. Cambridge University Press. Cambridge, UK. 
%% ~put references to the literature/web site here ~
}
\author{Long Qu 
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dt.int2}}, \code{\link{dt.sad}}, \code{\link{dt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% __ONLY ONE__ keyword per line
