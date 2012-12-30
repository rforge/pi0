get.p.dd.dw <-
function(y, permutedTrt, r=seq_len(ncol(y)), #test=F, 
                        #distObj=dist(y), 
                        #cperm.mat, 
                        wtmethod=0, 
                        eps=1e-8
                    )      
## b=permutation index; r=dimension index; 
{
    B=nperms.permutedTrt(permutedTrt)
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
#    b=if(test) 1:B else 1
#    ans=matrix(NA_real_, length(r), length(b))
    ans=numeric(length(r))
    N=nrow(y)
    #if(missing(cperm.mat)) cperm.mat=apply(permutedTrt,2,function(kk)(1:N)[-kk])

#    contrast.mat=matrix(0,choose(N,2),N); k=1
#    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
#    diffy=contrast.mat%*%y; diffy2=diffy*diffy
#    all.ddelta.dw=diffy2/pmax(1e-8,sqrt(rowSums(diffy2)))/2
##    all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,dist(y))/2  
    all.ddelta.dw=apply(y[,r,drop=FALSE],2L,dist)^2/dist(y)*0.5   ## these 2 lines replace the above 5 lines
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0

    for(r.i in seq(along=r)){
        dz.dw=.Call(mrppstats,all.ddelta.dw[,r.i],permutedTrt, as.numeric(wtmethod[1]), PACKAGE='MRPP')
        ans[r.i]=mean(dz.dw[1]-dz.dw>= -eps)
    }
    drop(ans)
}

