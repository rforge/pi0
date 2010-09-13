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
	    alpha.in, alpha.del = 0, stepwise = FALSE, 
	    size.in=0L, cperm.mat,...)
bsmrpp(y,perm.mat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, alpha.del=0, stepwise=FALSE, 
	 size.in=1L, cperm.mat, Bperm=ncol(perm.mat), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	Numeric data matrix, with columns being variables and rows being observations.
}
  \item{perm.mat}{
	A permutation of indices of the first treatment group. Each column represent a permutation.
}
  \item{verbose}{
	Logical or numeric scalar. Print messages every \code{verbose} iteration. \code{TRUE} is treated as 1. If \code{FALSE}, no messages are printed.
}
  \item{niter}{
	Maximum number of iterations. 
}
  \item{importance}{
	Either \code{'dp.dw'} or \code{'p.dd.dw'}. The importance measure. \code{'dp.dw'} means the derivative of approximate p-value with respect to weights. \code{'p.dd.dw'} means the permutation p-value using the derivative of Euclidean distance with respect to weights as the new distance measure. 
}
  \item{alpha.in}{
	Desired minimally allowable level of importance for selected variables. 
}
  \item{alpha.del}{
	Desired minimally allowable MRPP p-value for excluded variables. 
}
  \item{stepwise}{
	Logical. If \code{TRUE}, then only the variables with least importance are removed. Otherwise, all variables with importance less than \code{alpha.in} are removed.
}
  \item{size.in}{
	Minimally allowable number of selected variables. 
}
  \item{cperm.mat}{
	Similar to \code{perm.mat}, except that this is the permutation indices matrix for the second treatment group. 
}
  \item{Bperm}{
	Number of outer permutations. 
}
  \item{\dots}{
	Additional arguments passed to \code{mrpp.test.dist}, \code{get.dp.dw.kde} or \code{get.p.dd.dw}.
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
\code{back.search} returns a list, each element corresponding to results in each iteration. Each element is a list by itself, with elements
\item{iter}{Iteration index;}
\item{var.idx}{Variable index, in the order of variable importance;}
\item{influence}{Sorted ariable importance;}
\item{p.value}{MRPP p-value for selected variables;}
\item{deleted.p.value}{MRPP p-value for excluded variables.}

\code{bsmrpp} returns the last element of the result from \code{back.search}, except that
\item{raw.p.value}{replaces the original p-value from \code{back.search};}
\item{adj.p.value}{the actual p-value (adjusted) computed by nesting variable selection in permutations;}
\item{permuted.p.values}{a vector of raw.p.values, one for each outer permutation.}
}
\references{
Long Qu, Dan Nettleton, and Jack C. M. Dekkers. RELATIVE VARIABLE IMPORTANCE AND BACKWARD VARIABLE SELECTION FOR THE MULTIRESPONSE
PERMUTATION PROCEDURE, WITH APPLICATIONS TO HIGH DIMENSIONAL GENOMIC DATA. (submitted)
}
\author{
Long Qu
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
