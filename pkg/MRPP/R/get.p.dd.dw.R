get.p.dd.dw <-
function(y, permmat, r=1:ncol(as.matrix(y)), #test=F, 
                        distObj=dist(y), 
                        cpermmat, wtmethod=integer(1)
                    )      
## b=permutation index; r=dimension index; 
{
    B=ncol(permmat)
#    b=if(test) 1:B else 1
#    ans=matrix(NA_real_, length(r), length(b))
    ans=numeric(length(r))
    N=as.integer(nrow(y))
    if(missing(cpermmat)) cpermmat=apply(permmat,2,function(kk)(1:N)[-kk])

    contrast.mat=matrix(0,choose(N,2),N); k=1
    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
    all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2

    for(r.i in seq(along=r)){
        dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permmat,cpermmat,nrow(permmat),B,N,as.integer(wtmethod[1]),ans=double(B))$ans
        ans[r.i]=mean(dz.dw[1]-dz.dw>=-1e-8)
    }
    drop(ans)
}

