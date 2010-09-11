back.search <-
function(y,perm.mat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, stepwise=FALSE, cperm.mat, ...)
## y is a data matrix, with col's being variables and rows being observations
{
    if(!is.matrix(y))y=as.matrix(y)
    importance=match.arg(importance)
    if(missing(alpha.in)) alpha.in=if(importance=='dp.dw') 0 else 0.1
    N=nrow(y)
    if(missing(cperm.mat)) cperm.mat=apply(perm.mat,2,function(kk)(1:N)[-kk])
    B=ncol(perm.mat)

    ans=vector('list')
    idx=1:ncol(y)
    i=1
    repeat{
        if(verbose) {cat('iteration',i-1,'...')
                    time0=proc.time()[3]}
        dist0=dist(y[,idx,drop=FALSE])
        mrpp.stats0=mrpp.test.dist(dist0,perm.mat=perm.mat,...)$all.stat
        imptnc=if(importance=='dp.dw') 
                    get.dp.dw.kde(y[,idx,drop=FALSE],perm.mat,distObj=dist0,mrpp.stats=mrpp.stats0, cperm.mat=cperm.mat,...) 
               else get.p.dd.dw(y[,idx,drop=FALSE],perm.mat,distObj=dist0,cperm.mat=cperm.mat,...)
        var.rank=order(imptnc)
        ans[[i]]=list(iter=i-1, var.idx=idx[var.rank], influence=imptnc[var.rank],
                      p.value=mean(mrpp.stats0[1]>=mrpp.stats0-min(c(1e-8,.5/B))),
                      deleted.p.value=NA_real_)
        if(alpha.del>0) {
            dist.del=dist(y[,-idx,drop=FALSE])
            if(all(!is.na(dist.del)))
                ans[[i]]$deleted.p.value=mrpp.test.dist(dist.del, perm.mat=perm.mat)$p.value
        }
        if(verbose) {
          cat('\b\b\b:\t',length(idx),'genes left; mrpp.p =',ans[[i]]$p.value,';', 
                        'deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3]-time0,'seconds passed;',fill=TRUE)
        }
        if(all(imptnc<alpha.in) || i-1>=niter || isTRUE(ans[[i]]$deleted.p.value<=alpha.del)) return(ans)
        i=i+1
        if(stepwise) idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.in]
        if(length(idx)==0) {
            warning('not converged'); 
            ans[[i]]=list(iter=i-1, var.idx=numeric(0), influence=numeric(0),
                      p.value=numeric(0),
                      deleted.p.value=ans[[1]]$p.value)
            return(ans)}
    }
}

