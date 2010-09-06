back.search <-
function(y,permmat, verbose=TRUE, niter=Inf, 
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
        mrpp.stats0=apply(permmat,2,get.mrpp.stat,dist.mat=dist0)
        imptnc=if(importance=='dp.dw') get.dp.dw.kde(y[,idx,drop=FALSE],permmat,dist.mat=dist0,mrpp.stats=mrpp.stats0, ...) 
               else get.p.dd.dw(y[,idx,drop=FALSE],permmat,dist.mat=dist0,...)
        ans[[i]]=list(iter=i-1, var.idx=idx[order(imptnc)], influence=sort(imptnc),
                      p.value=mean(mrpp.stats0[1]>=mrpp.stats0),
                      deleted.p.value=NA_real_)
        if(verbose) {
          if(alpha.del>0) {
            dist.del=as.matrix(dist(y[,-idx,drop=FALSE]))
            mrpp.stats.del=apply(permmat,2,get.mrpp.stat,dist.mat=dist.del)
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

