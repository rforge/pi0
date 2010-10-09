back.search <-
function(y,perm.mat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, stepwise=FALSE, size.in=0L, cperm.mat, ...)
## y is a data matrix, with col's being variables and rows being observations
{
    if(!is.matrix(y))y=as.matrix(y)
    importance=match.arg(importance)
    if(missing(alpha.in)) alpha.in=if(importance=='dp.dw') 0 else 0.1
    N=nrow(y)
    if(missing(cperm.mat)) cperm.mat=apply(perm.mat,2,function(kk)(1:N)[-kk])
    B=ncol(perm.mat)

    ans=vector('list')
    R=ncol(y)
    idx=1:R  # inclusion set
    xcl=integer(0)
    i=1
    imptnc.threshold=Inf
    imptnc=rep(-Inf, R)
    next.deleted.p=NA_real_
    repeat{
        if(verbose && isTRUE(i%%verbose==0)) {cat('iteration',i-1,'...')
                    time0=proc.time()[3]}
        idx=idx[imptnc<imptnc.threshold]
        if(length(idx)==0){
            ans[[i]]=list(iter=i-1,var.idx=integer(0), influence=numeric(0), 
                        p.value=numeric(0),deleted.p.value=ans[[1]]$p.value)
            if(verbose)message('no variables left in the inclusion set')
            return(ans)
        }
        dist0=dist(y[,idx,drop=FALSE])
        mrpp.stats0=mrpp.test.dist(dist0,perm.mat=perm.mat,cperm.mat=cperm.mat,...)$all.stat
        imptnc=if(importance=='dp.dw') 
                    get.dp.dw.kde(y[,idx,drop=FALSE],perm.mat,distObj=dist0,mrpp.stats=mrpp.stats0, cperm.mat=cperm.mat,...) 
               else get.p.dd.dw(y[,idx,drop=FALSE],perm.mat,cperm.mat=cperm.mat,...)
        var.rank=order(imptnc)
        ans[[i]]=list(iter=i-1, var.idx=idx[var.rank], influence=imptnc[var.rank],
                      p.value=mean(mrpp.stats0[1]>=mrpp.stats0-min(c(1e-8,.5/B))),
                      deleted.p.value=next.deleted.p)
        imptnc.threshold=if(stepwise) max(imptnc) else alpha.in
#        idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.in]
        if(alpha.del>=0) {
            xcl=c(xcl, idx[imptnc>=imptnc.threshold])
            dist.del=dist(y[,xcl,drop=FALSE])
            if(all(!is.na(dist.del)) && length(xcl)>0)
#                ans[[i]]$deleted.p.value=mrpp.test.dist(dist.del, perm.mat=perm.mat)$p.value
                next.deleted.p=mrpp.test.dist(dist.del, perm.mat=perm.mat,cperm.mat=cperm.mat)$p.value
        }
        if(verbose && isTRUE(i%%verbose==0)) {
          cat('\b\b\b:\t',length(idx),'genes left; mrpp.p =',ans[[i]]$p.value,';', 
                        'deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3]-time0,'seconds passed;',fill=TRUE)
        }
        if( all(imptnc<=alpha.in) || i-1>=niter || 
#            isTRUE(ans[[i]]$deleted.p.value<=alpha.del) ||
            isTRUE(next.deleted.p<=alpha.del) ||
            isTRUE(R-length(xcl)<size.in) ) return(ans)
        i=i+1
#        if(stepwise) idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.in]
#        if(length(idx)==0) {
#            warning('not converged'); 
#            ans[[i]]=list(iter=i-1, var.idx=numeric(0), influence=numeric(0),
#                      p.value=numeric(0),
#                      deleted.p.value=ans[[1]]$p.value)
#            return(ans)}
    }
}

