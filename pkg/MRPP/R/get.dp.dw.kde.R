get.dp.dw.kde <-
function(y, perm.mat, r=1:ncol(as.matrix(y)), test=FALSE, 
        distObj=dist(y), 
        mrpp.stats=mrpp.test.dist(distObj,perm.mat=perm.mat,wtmethod=wtmethod[1])$all.statistics,
        bw=bw.mse.pdf.asym(mrpp.stats,iter.max=1L), #cperm.mat, 
        wtmethod=0)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
    min.wts=0   ### CHECK ME: I cannot remember why this weights were introduced. Set it to zero for now see what problems show up...
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    ans=matrix(NA_real_, length(b), length(r))
    if(!is.matrix(y))y=as.matrix(y)
    N=as.integer(nrow(y))
    #if(missing(cperm.mat)) cperm.mat=apply(perm.mat,2,function(kk)(1:N)[-kk])

#    weight=matrix(NA_real_, B, length(b))   ## this may require large memory when test=TRUE
#    for(b.i in 1:length(b))
#      weight[,b.i]=pmax(min.wts,dnorm((mrpp.stats[b[b.i]]-mrpp.stats),0,bw))
    weight = dnorm(outer(mrpp.stats, mrpp.stats[b], '-'), 0, bw)    ## this replaces the above 3 lines


#    contrast.mat=matrix(0,choose(N,2),N); k=1
#    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
#    #all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero
#    all.ddelta.dw=(contrast.mat%*%y)^2/distObj/2 ## when denom is zero, the numerator is also zero. 
#        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
    all.ddelta.dw=apply(y,2L,dist)^2/distObj*0.5   ## these 2 lines replace the above 5 lines
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0

    for(r.i in seq(along=r)){
        #dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],perm.mat,cperm.mat,nrow(perm.mat),B,N,as.integer(wtmethod[1]),ans=double(B),PACKAGE='MRPP',DUP=FALSE)$ans
        dz.dw=.Call('mrppstats', all.ddelta.dw[,r.i], perm.mat, as.numeric(wtmethod), PACKAGE='MRPP')
#        for(b.i in 1:length(b)){
#            dd.dw=dz.dw[b[b.i]]-dz.dw
#            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/B  #length(b)
#        }
        ans[, r.i] = colMeans(weight * outer(-dz.dw, dz.dw[b], '+'))    ## this lines replace the above 3 lines
    }
    drop(ans)
}

