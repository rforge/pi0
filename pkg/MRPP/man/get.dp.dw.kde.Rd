\name{get.dp.dw.kde}
\Rdversion{1.1}
\alias{get.dp.dw.kde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute variable importance measure of partial derivative of p-values w.r.t. weights
}
\description{
Compute variable importance measure of partial derivative of p-values w.r.t. weights
}
\usage{
get.dp.dw.kde(y, permmat, r = 1:ncol(as.matrix(y)), test = FALSE, distObj = dist(y), mrpp.stats = mrpp.test.dist(distObj, perm.mat = permmat,wtmethod=wtmethod[1])$all.stat, bw = bw.mse.pdf.asym(mrpp.stats), cpermmat, wtmethod=integer(1), min.wts = 1e-08)
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
  \item{test}{
%%     ~~Describe \code{test} here~~
}
  \item{distObj}{
%%     ~~Describe \code{distObj} here~~
}
  \item{mrpp.stats}{
%%     ~~Describe \code{mrpp.stats} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{cpermmat}{
%%     ~~Describe \code{cpermmat} here~~
}
  \item{wtmethod}{
%%     ~~Describe \code{cpermmat} here~~
}
  \item{min.wts}{
%%     ~~Describe \code{min.wts} here~~
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
function(y, permmat, r=1:ncol(as.matrix(y)), test=F, 
                        distObj=dist(y), 
                        mrpp.stats=mrpp.test.dist(distObj,perm.mat=permmat)$all.stat,
                        bw=bw.mse.f.asym(mrpp.stats), cpermmat,
                        min.wts=1e-8)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
    B=length(mrpp.stats)
    b=if(test) 1:B else 1
    ans=matrix(NA_real_, length(r), length(b))
    if(!is.matrix(y))y=as.matrix(y)
    N=as.integer(nrow(y))
    if(missing(cpermmat)) cpermmat=apply(permmat,2,function(kk)(1:N)[-kk])

    weight=matrix(NA_real_, B, length(b))
    for(b.i in 1:length(b))
      weight[,b.i]=pmax(min.wts,dnorm((mrpp.stats[b[b.i]]-mrpp.stats),0,bw))

    contrast.mat=matrix(0,choose(N,2),N); k=1
    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
    all.ddelta.dw=abs(contrast.mat\%*\%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero

    for(r.i in seq(along=r)){
        dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permmat,cpermmat,nrow(permmat),B,N,ans=double(B))$ans
        for(b.i in 1:length(b)){
            dd.dw=dz.dw[b[b.i]]-dz.dw
            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/length(b)
        }
    }
    drop(ans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
