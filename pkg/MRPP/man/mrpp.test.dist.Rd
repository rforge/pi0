\name{mrpp.test.dist}
\Rdversion{1.1}
\alias{mrpp.test.dist}
\alias{mrpp.test.matrix}
\alias{mrpp.test.formula}
\alias{mrpp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mrpp.test.dist(y, trt, B = choose(length(trt), table(trt)[1]), perm.mat, wtmethod=0)
mrpp.test.formula(y, B, perm.mat)
mrpp.test.matrix(y, ...)
mrpp.test(y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{trt}{
%%     ~~Describe \code{trt} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{perm.mat}{
%%     ~~Describe \code{perm.mat} here~~
}
  \item{wtmethod}{
%%     ~~Describe \code{perm.mat} here~~
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
function(y, trt, B=choose(length(trt),table(trt)[1]), perm.mat) ## this is C code
## y is a dist object
{
    if(missing(y) || !inherits(y,'dist')) stop('dist object missing or incorrect')
    N=attr(y,'Size')
    if(missing(perm.mat)) {
        perm.mat=get.perm.mat(trt,B)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permutation matrix', deparse(substitute(perm.mat)))
    B=ncol(perm.mat)
    cperm.mat=apply(perm.mat, 2, function(kk)(1:N)[-kk])
    stats=.C('mrppstats',y,perm.mat,cperm.mat,nrow(perm.mat),B,N,ans=double(B))$ans
    ans=list(statistic=c("MRPP statistic"=stats[1]), all.statistics=stats, 
             p.value=mean(stats[1]>=stats), parameter=c("number of permutations"=B),
             data.name=dname, .Random.seed=attr(perm.mat,'.Random.seed'),
             method='2-sample MRPP test')
    class(ans)='htest'
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
