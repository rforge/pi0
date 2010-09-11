\name{back.search}
\Rdversion{1.1}
\alias{back.search}
\alias{bsmrpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backward variable selection for MRPP
}
\description{
Backward variable selection for MRPP
}
\usage{
back.search(y, perm.mat, verbose = TRUE, niter = Inf, 
	    importance = c("dp.dw", "p.dd.dw"), 
	    alpha.in, alpha.del = 0, stepwise = FALSE, cperm.mat,...)
bsmrpp(y,perm.mat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, alpha.del=0, stepwise=FALSE, cperm.mat, Bperm=ncol(perm.mat), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{perm.mat}{
%%     ~~Describe \code{perm.mat} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{importance}{
%%     ~~Describe \code{importance} here~~
}
  \item{alpha.in}{
%%     ~~Describe \code{alpha.in} here~~
}
  \item{alpha.del}{
%%     ~~Describe \code{alpha.del} here~~
}
  \item{stepwise}{
%%     ~~Describe \code{stepwise} here~~
}
  \item{cperm.mat}{
%%     ~~Describe \code{stepwise} here~~
}
  \item{Bperm}{
%%     ~~Describe \code{stepwise} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function(y,perm.mat, verbose=TRUE, niter=Inf, 
                     importance=c('dp.dw','p.dd.dw'),
                     alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
                     alpha.del=0, stepwise=FALSE, ...)
## y is a data matrix, with col's being variables and rows being observations
{
    ans=vector('list')
    importance=match.arg(importance)
    if(missing(alpha.in)) alpha.in=if(importance=='dp.dw') 0 else 0.1
    idx=1:ncol(y)
    i=1
    repeat{
        if(verbose) cat('iteration',i-1,'...')
        time0=proc.time()[3]
        dist0=as.matrix(dist(y[,idx,drop=FALSE]))
        mrpp.stats0=apply(perm.mat,2,get.mrpp.stat,dist.mat=dist0)
        imptnc=if(importance=='dp.dw') get.dp.dw.kde(y[,idx,drop=FALSE],perm.mat,dist.mat=dist0,mrpp.stats=mrpp.stats0, ...) 
               else get.p.dd.dw(y[,idx,drop=FALSE],perm.mat,dist.mat=dist0,...)
        ans[[i]]=list(iter=i-1, var.idx=idx[order(imptnc)], influence=sort(imptnc),
                      p.value=mean(mrpp.stats0[1]>=mrpp.stats0),
                      deleted.p.value=NA_real_)
        if(verbose) {
          if(alpha.del>0) {
            dist.del=as.matrix(dist(y[,-idx,drop=FALSE]))
            mrpp.stats.del=apply(perm.mat,2,get.mrpp.stat,dist.mat=dist.del)
            ans[[i]]$deleted.p.value=mean(mrpp.stats.del[1]>=mrpp.stats.del)
          }
          cat('\b\b\b:\t',length(idx),'genes left; mrpp.p =',ans[[i]]$p.value,';', 
                        'deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3]-time0,'seconds passed;',fill=TRUE)
        }
        if(all(imptnc<alpha.in) || i-1>=niter || isTRUE(ans[[i]]$deleted.p.value<alpha.del)) return(ans)
        i=i+1
        if(stepwise) idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.in]
        if(length(idx)==0) {warning('not converged'); return(ans)}
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
