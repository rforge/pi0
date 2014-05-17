\name{dlassot}
\alias{dlassot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density function for the lassot prior distribution
}
\description{
Density function for the lassot prior distribution
}
\usage{
dlassot(x, rate = 1/continuityFactor/ffp1df, alpha = 2/ffp1df, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector of quantiles.
}
  \item{rate}{Numeric vector of rates (lambda).
}
  \item{alpha}{Numeric vector of alpha values.
}
  \item{log}{Logical; if \code{TRUE}, probabilities p are given as log(p).
}
}
\details{
Lassot distribution has a central similar to double-exponential distribution and tails similar to scaled central t distribution.
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
\author{Long Qu
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
\keyword{ distribution }
