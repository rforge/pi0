\name{dt.sad and pt.sad}
\Rdversion{1.1}
\alias{pt.sad}
\alias{dt.sad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Saddle Point Approximation of noncentral t-distribution
}
\description{ Density and cumulative distribution function of noncentral t-distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.sad(x, df, ncp = 0, log = FALSE, normalize = c("approximate", "derivative", "integral", "none"), epsilon = 1e-04)
pt.sad(q, df, ncp = 0, log = FALSE, epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{ numeric vector of quantiles
%%     ~~Describe \code{q} here~~
}
  \item{df}{ numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{ numeric vector of noncentrality parameter
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{ logical; whether \code{log} should be taken. 
%%     ~~Describe \code{log} here~~
}
  \item{epsilon}{ a small numeric scalar; if the difference between \code{q} and \code{ncp} is closer than this, results will be computed differently. 
%%     ~~Describe \code{epsilon} here~~
}
  \item{normalize}{ the way to normalize the approximate density so that it is closer to a true density. 
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{ \code{dt.sad} returns density; \code{pt.sad} returns the probability. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Broda, Simon and Paolella, Marc S. 2007. Saddlepoint approximations for the doubly noncentral t distribution, Computational Statistics & Data Analysis, 51,6, 2907-2918.
%% ~put references to the literature/web site here ~
}
\author{Long Qu 
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{dt.int2}}, \code{\link{dt.lap}}, \code{\link{dt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
%\keyword{ univar }% __ONLY ONE__ keyword per line
