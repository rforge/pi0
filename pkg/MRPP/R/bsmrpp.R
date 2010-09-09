bsmrpp<-
function(y,permmat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, stepwise=FALSE, cpermmat, Bperm=ncol(permmat), ...)
{   if(!is.matrix(y)) y=as.matrix(y)
    N=nrow(y)
    if(missing(cpermmat)) cpermmat=apply(permmat,2,function(kk)(1:N)[-kk])
    selected.pvals=numeric(Bperm)
    bsfit=tail(back.search(y,permmat, FALSE, niter, importance, alpha.in, alpha.del, stepwise, cpermmat,...),1)[[1]]
    selected.pvals[1]=bsfit$p.value
    lastperm=perm1=permmat[,1]; lastcperm=cperm1=cpermmat[,1]
    for(b in 2:Bperm) {
        if(verbose) cat("iteration:",b-1," out of",Bperm,"\t\t\r")

        permmat[,b-1]=lastperm; cpermmat[,b-1]=lastcperm
        lastperm=permmat[,1]=permmat[,b]; lastcperm=cpermmat[,1]=cpermmat[,b]
        permmat[,b]=perm1;  cpermmat[,b]=cperm1

        selected.pvals[b]=tail(back.search(y,permmat, FALSE, niter, importance, 
                          alpha.in, alpha.del, stepwise, cpermmat,...),1)[[1]]$p.value
    }
    bsfit$raw.p.value=selected.pvals[1]
    bsfit$p.value=NULL
    bsfit$adj.p.value=mean(selected.pvals[1]-selected.pvals>=-min(c(1e-8,.5/Bperm)))
    bsfit
}
