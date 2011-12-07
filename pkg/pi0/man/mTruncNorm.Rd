\name{mTruncNorm}
\Rdversion{1.1}
\alias{mTruncNorm}
\alias{mTruncNorm.int2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moments of truncated normal distribution and the integral in the noncentral t-distribution
}
\description{Compute the moments of truncated normal distribution and the integral that appears in the noncentral t-distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mTruncNorm(r = 1, mu = 0, sd = 1, lower = -Inf, upper = Inf, 
        approximation = c("int2", "laplace", "numerical"), 
        integral.only = FALSE, ...)
mTruncNorm.int2(r = as.integer(1), mu = 0, sd = 1, lower = -Inf, 
        upper = Inf, takeLog = TRUE, ndiv = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{the order of moments to be computed. It could be noninteger, but has to be nonnegative. This is also the degrees of freedom for the noncentral t-distribution.
%%     ~~Describe \code{r} here~~
}
  \item{mu}{mean of the normal distribution, before truncating. 
%%     ~~Describe \code{mu} here~~
}
  \item{sd}{SD of the normal distribution, before truncating. 
%%     ~~Describe \code{sd} here~~
}
  \item{lower}{lower truncation point
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{upper truncation point
%%     ~~Describe \code{upper} here~~
}
  \item{approximation}{Method of approximation. \code{int2} is exact for \emph{int}eger \code{r} and \emph{int}erpolate to noninteger \code{r}. 
\code{laplace} uses laplacian approximation. \code{numerical} uses nuemerical integration. 
%%     ~~Describe \code{approximation} here~~
}
\item{integral.only}{logical. If \code{TRUE}, only the integral in noncentral t-distribution is returned. Otherwise, it is normalized to be the rth moments of truncated normal distribution. 
%%     ~~Describe \code{integral.only} here~~
}
  \item{takeLog}{logical. If \code{TRUE} and \code{r} is not an integer, the polyomial interpolation will be on the log scale. But final result is on the original scale.
}
  \item{ndiv}{number of points with closes integer \code{r} to be used in polynomial interpolation. 
}
  \item{\dots}{other arguments passed to \code{mTruncNorm.int2}
%%     ~~Describe \code{\dots} here~~
}
}
\details{\code{mTruncNorm.int2} uses iterative relation over \code{r} to compute the integral iteratively starting from \code{r=0} and \code{r=1} whose analytic results are available. 
If \code{r} is not an integer, the nearest \code{ndiv} nonnegative integer \code{r} will be used to do divided difference polynomial interpolation. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{numeric vector. If \code{integral.only} is \code{TRUE}, this is the integral in the noncentral t-density; otherwise this is the rth moments of truncated normal distribution. 
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

\seealso{\code{\link{dt}}, \code{\link{pt}}, \code{\link{dt.int2}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ distribution }% __ONLY ONE__ keyword per line
