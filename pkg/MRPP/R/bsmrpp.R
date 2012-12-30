mrppBVS.test<-
function(y,permutedTrt, Bperm=nperms.permutedTrt(permutedTrt), 
         importance=c('dp.dw','p.dd.dw'),
         alpha.in, alpha.del=0, 
	 size.in=1L, stepwise=FALSE, verbose=TRUE, niter=Inf, ...)
{
    dname=paste("Response data", deparse(substitute(y)), 'and permuted treatment', deparse(substitute(permutedTrt)))
    if(!is.matrix(y)) y=as.matrix(y)
    N=nrow(y)
    B=nperms.permutedTrt(permutedTrt)
    if(Bperm > B) stop("Current implementation requires Bperm <= nperms.permutedTrt(permutedTrt)");
    ntrt=ntrt.permutedTrt(permutedTrt)
    importance=match.arg(importance)

#    if(missing(cperm.mat)) cperm.mat=apply(permutedTrt,2,function(kk)(1:N)[-kk])
    if(missing(alpha.in)) alpha.in=if(importance=='dp.dw') 0 else 0.1
    selected.pvals=numeric(Bperm)
    fit0=mrppBVS(y=y,permutedTrt=permutedTrt, verbose=FALSE, niter=niter, importance=importance, 
            alpha.in=alpha.in, alpha.del=alpha.del, stepwise=stepwise, size.in=size.in, 
            ...)
    bsfit=tail(fit0,1L)[[1L]]
    selected.pvals[1]=if(length(bsfit$p.value)>0) bsfit$p.value else 1
#    lastperm=perm1=permutedTrt[,1]; lastcperm=cperm1=cperm.mat[,1]

    if(Bperm < B) {
        bseqs=sample(B-1L, Bperm-1L)+1L
    }else if (B == 1L) {
        bseqs=integer(0L)
    }else bseqs=2:B
    for(b.i in seq_along(bseqs)) {
        if(verbose && isTRUE(b.i%%verbose==0L)) cat("outer permutation:",b.i-1L," out of",Bperm,"\t\t\r")

#        permutedTrt[,b-1L]=lastperm; cperm.mat[,b-1]=lastcperm
#        lastperm=permutedTrt[,1]=permutedTrt[,b]; lastcperm=cperm.mat[,1]=cperm.mat[,b]
#        permutedTrt[,b]=perm1;  cperm.mat[,b]=cperm1
        if(is.na(attr(permutedTrt, 'idx')[1L])) {
            for(tt in seq(ntrt)) {
                permutedTrt[[tt]][,1L] -> tmp
                permutedTrt[[tt]][,1L] <- permutedTrt[[tt]][,bseqs[b.i]]
                tmp                    -> permutedTrt[[tt]][,bseqs[b.i]]
            }
        }else {
                attr(permutedTrt, 'idx')[1L] -> tmp
                attr(permutedTrt, 'idx')[1L] <- attr(permutedTrt, 'idx')[bseqs[b.i]]
                tmp                          -> attr(permutedTrt, 'idx')[bseqs[b.i]]
        }

        selected.pvals[b.i+1L]={tmp=tail(mrppBVS(y=y,permutedTrt=permutedTrt, verbose=FALSE, niter=niter, 
                        importance=importance, alpha.in=alpha.in, alpha.del=alpha.del, 
                        stepwise=stepwise, size.in=size.in, ...),1L)[[1L]]$p.value;
                           if(length(tmp)>0) tmp else 1}
    }
    bsfit$statistic=c('MRPP backward selected raw p-value'=selected.pvals[1])
    bsfit$p.value=mean(selected.pvals - selected.pvals[1] < 0.5 / B)
    bsfit$all.statistics=selected.pvals
    bsfit$method="Permutation test with MRPP backward variable selection"
    parms=attr(fit0, 'parameter')
    bsfit$parameter=c(importance=match(parms$importance, c('dp.dw','p.dd.dw')), unlist(parms[-1]), 
        '#selected variables'=length(bsfit$var.idx))
    bsfit$data.name=dname
    class(bsfit)='htest'
    bsfit
}
