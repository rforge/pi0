\name{dlassot}
\alias{dlassot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dlassot(x, rate = 1/continuityFactor/ffp1df, alpha = 2/ffp1df, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{rate}{
%%     ~~Describe \code{rate} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (x, rate = 1/continuityFactor/ffp1df, alpha = 2/ffp1df, 
    log = FALSE) 
{
    a = alpha
    stopifnot(all(rate > 0, alpha > 1/ffp1df))
    lambda = rate
    f = lassotFactor
    const = 1/((2 - 2 * exp(-a * f * lambda))/lambda + 1/(-f + 
        a * (1 + f^2) * lambda) * 2 * a * f^(2 - (a * (1 + f^2) * 
        lambda)/f) * hypergeo::hypergeo((a * (1 + f^2) * lambda)/(2 * 
        f), 1/2 * (-1 + (a * lambda)/f + a * f * lambda), (f + 
        a * (1 + f^2) * lambda)/(2 * f), -(1/f^2)))
    logans = log(Re(const)) + ifelse(abs(x) <= lassotFactor * 
        a, -abs(x) * rate, -ffp1df/2 * rate * a * log(1 + x * 
        x/a/a))
    if (isTRUE(log)) 
        logans
    else exp(logans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
