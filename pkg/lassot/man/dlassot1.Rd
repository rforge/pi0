\name{dlassot1}
\alias{dlassot1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dlassot1(x, a, m, df = 2 * m - 1, s = a/sqrt(df), log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, a, m, df = 2 * m - 1, s = a/sqrt(df), log = FALSE) 
{
    half.const = pt(-lassotFactor * sqrt(df), df) * a * beta(df/2, 
        0.5) + a * lassotFactor/m/log1pff * (1 - (1 + lassotFactor * 
        lassotFactor)^-m)
    logans = ifelse(abs(x) <= lassotFactor * a, -abs(x) * m * 
        log1pff/lassotFactor/a, -m * log(1 + x * x/a/a)) - log(half.const) - 
        log(2)
    if (isTRUE(log)) 
        logans
    else exp(logans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
