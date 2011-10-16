\name{dt.int2}
\Rdversion{1.1}
\alias{dt.int2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{INTerpolation of INTeger degrees of freedom noncentral t-density
}
\description{This funciton evaluates the noncentral t-density using an iterative procedure for integer degrees of freedom. This is much faster than two calls to the \code{pt} approach. 
For non-integer degrees of freedom, the polynomial interpolation is used to approximate the density. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dt.int2(x, df, ncp, log = FALSE, ndiv = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector of quantiles
%%     ~~Describe \code{x} here~~
}
  \item{df}{A numeric vector degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{ncp}{A numeric vector of noncentrality parameters
%%     ~~Describe \code{ncp} here~~
}
  \item{log}{logical; if \code{TRUE}, log densities are returned. 
%%     ~~Describe \code{log} here~~
}
  \item{ndiv}{numeric; the number of points used for polynomial interpolation
%%     ~~Describe \code{ndiv} here~~
}
}
\details{This function uses the iterative relation for the integral in the noncentral t-density. It starts with \code{df=0} and \code{df=1}, 
and then iteratively computes the integral for larger \code{df}. For non-integer \code{df}, it uses \code{ndiv} nearest points to perform a 
divided difference polynomial interpolation approximation. For integer \code{df}, this function is about 2 to 3 times faster than \code{dt} function and is exact.  
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric vector of densities. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%Qu, Nettleton, Dekkers. Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. (unpublished)
%}
\author{Long Qu \email{lqu@biostatsolutions.com}
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dtn.mix}}, \code{\link{dt}}, \code{\link{dt.lap}}, \code{\link{dt.sad}}, \code{\link{mTruncNorm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% __ONLY ONE__ keyword per line
