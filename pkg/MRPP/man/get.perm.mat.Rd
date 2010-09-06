\name{get.perm.mat}
\Rdversion{1.1}
\alias{get.perm.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.perm.mat(trt, B = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{
%%     ~~Describe \code{trt} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function(trt, B=100) ## trt needs have exactly 2 levels; other situations are not implemented yet
{
    if(!exists('.Random.seed', envir=globalenv())) runif(1)
    save.seed=get('.Random.seed', envir=globalenv())
    n=table(trt)
    perms=replicate(B, sample(sum(n),min(n)))
    perms[,1]=which(as.factor(trt)==names(n[which.min(n)]))
    attr(perms,'.Random.seed')=save.seed
    perms
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
