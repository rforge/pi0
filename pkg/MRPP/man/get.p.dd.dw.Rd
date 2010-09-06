\name{get.p.dd.dw}
\Rdversion{1.1}
\alias{get.p.dd.dw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute variable importance measure by using partial derivative of distance function w.r.t. weights as the new distance function
}
\description{
Compute variable importance measure by using partial derivative of distance function w.r.t. weights as the new distance function
}
\usage{
get.p.dd.dw(y, permmat, r = 1:ncol(as.matrix(y)), distObj = dist(y), cpermmat, wtmethod=integer(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{permmat}{
%%     ~~Describe \code{permmat} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{distObj}{
%%     ~~Describe \code{distObj} here~~
}
  \item{cpermmat}{
%%     ~~Describe \code{cpermmat} here~~
}
  \item{wtmethod}{
%%     ~~Describe \code{cpermmat} here~~
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
function(y, permmat, r=1:ncol(as.matrix(y)), #test=F, 
                        distObj=dist(y), 
                        cpermmat
                    )      
## b=permutation index; r=dimension index; 
{
    B=ncol(permmat)
#    b=if(test) 1:B else 1
#    ans=matrix(NA_real_, length(r), length(b))
    ans=numeric(length(r))
    N=as.integer(nrow(y))
    if(missing(cpermmat)) cpermmat=apply(permmat,2,function(kk)(1:N)[-kk])

    contrast.mat=matrix(0,choose(N,2),N); k=1
    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
    all.ddelta.dw=abs(contrast.mat\%*\%y)^2/pmax(1e-8,distObj)/2

    for(r.i in seq(along=r)){
        dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permmat,cpermmat,nrow(permmat),B,N,ans=double(B))$ans
        ans[r.i]=mean(dz.dw[1]-dz.dw>=-1e-8)
    }
    drop(ans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
