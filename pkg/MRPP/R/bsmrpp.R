bsmrpp<-
function(y,perm.mat, verbose=TRUE, niter=Inf, 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, #=if(match.arg(importance)=='dp.dw') 0 else 0.1, 
         alpha.del=0, stepwise=FALSE, size.in=1L, cperm.mat, Bperm=ncol(perm.mat), ...)
{   if(!is.matrix(y)) y=as.matrix(y)
    N=nrow(y)
    if(missing(cperm.mat)) cperm.mat=apply(perm.mat,2,function(kk)(1:N)[-kk])
    if(missing(alpha.in)) alpha.in=if(importance=='dp.dw') 0 else 0.1
    selected.pvals=numeric(Bperm)
    bsfit=tail(back.search(y=y,perm.mat=perm.mat, verbose=FALSE, niter=niter, importance=importance, 
            alpha.in=alpha.in, alpha.del=alpha.del, stepwise=stepwise, size.in=size.in, 
            cperm.mat=cperm.mat,...),1)[[1]]
    selected.pvals[1]=if(length(bsfit$p.value)>0) bsfit$p.value else 1
    lastperm=perm1=perm.mat[,1]; lastcperm=cperm1=cperm.mat[,1]
    for(b in 2:Bperm) {
        if(verbose) cat("iteration:",b-1," out of",Bperm,"\t\t\r")

        perm.mat[,b-1]=lastperm; cperm.mat[,b-1]=lastcperm
        lastperm=perm.mat[,1]=perm.mat[,b]; lastcperm=cperm.mat[,1]=cperm.mat[,b]
        perm.mat[,b]=perm1;  cperm.mat[,b]=cperm1

        selected.pvals[b]={tmp=tail(back.search(y=y,perm.mat=perm.mat, verbose=FALSE, niter=niter, 
                        importance=importance, alpha.in=alpha.in, alpha.del=alpha.del, 
                        stepwise=stepwise, size.in=size.in, cperm.mat=cperm.mat,...),1)[[1]]$p.value;
                           if(length(tmp)>0) tmp else 1}
    }
    bsfit$raw.p.value=selected.pvals[1]
    bsfit$p.value=NULL
    bsfit$adj.p.value=mean(selected.pvals[1]-selected.pvals>=-min(c(1e-8,.5/Bperm)))
    bsfit
}
