get.dp.dw.kde <-
function(y, permmat, r=1:ncol(as.matrix(y)), test=FALSE, 
                        distObj=dist(y), 
                        mrpp.stats=mrpp.test.dist(distObj,perm.mat=permmat,wtmethod=wtmethod[1])$all.stat,
                        bw=bw.mse.pdf.asym(mrpp.stats), cpermmat, wtmethod=integer(1),
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
    all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero

    for(r.i in seq(along=r)){
        dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permmat,cpermmat,nrow(permmat),B,N,as.integer(wtmethod[1]),ans=double(B),PACKAGE='MRPP')$ans
        for(b.i in 1:length(b)){
            dd.dw=dz.dw[b[b.i]]-dz.dw
            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/length(b)
        }
    }
    drop(ans)
}

