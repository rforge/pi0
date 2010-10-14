get.perm.mat <-
function(trt, B=100) ## trt needs have exactly 2 levels; other situations are not implemented yet
{
    if(!exists('.Random.seed', envir=globalenv())) runif(1)
    save.seed=get('.Random.seed', envir=globalenv())
    n=table(trt)
    perms=replicate(B, sort(sample(sum(n),min(n))))
    perms[,1]=which(as.factor(trt)==names(n[which.min(n)]))
    attr(perms,'.Random.seed')=save.seed
    perms
}


get.perm.mat <-
function(trt, B=100) ## permutation matrix for one way design
## returns a length(trt) by B matrix, if B is not larger than multinomial coefficient. 
{
    if(!exists('.Random.seed', envir=globalenv())) runif(1)
    save.seed=get('.Random.seed', envir=globalenv())
    n=table(trt)
    N=sum(n); 
#    multcoeff=exp(lfactorial(N)-sum(lfactorial(n)))
#    if(multcoeff>=B){   ## list all effective permutations
#        B=multcoeff
#        get.combn=function(vec, n){
#            if(length(vec)==n[1]) return(matrix(vec))
#            this.ans=combn(length(vec), n[1])
#            ans=apply(this.ans, 2, function(x){
#                    tmp=get.combn(vec[-x],n[-1])
#                    rbind(matrix(rep(vec[x],ncol(tmp)),n[1]), tmp)
#                })
#            matrix(ans,length(vec))
#        }
#        ans.mat=get.combn(N,n)
#    }else{  ## generate ramdom samples from all permutations
        get.choose=function(N,n)exp(lfactorial(N)-sum(lfactorial(n)))
        get.combn=function(vec, n, R){
            if(length(vec)==n[1]) return(matrix(vec))
            k=choose(length(vec), n[1])
            rest=get.choose(length(vec)-n[1],n[-1])
            if(k<=R) { 
                this.ans=combn(length(vec), n[1])
                if(rest*k>R){
                    sampled=table(sort(sample(rest*k, R)-1)%/%rest+1)
                    remainder=numeric(k)
                    remainder[as.numeric(names(sampled))]=sampled
                }else   remainder=rep(rest,k)
            }else {
                this.ans=combn2R(length(vec), n[1],R=R)
                remainder=rep(1L, R)
            }
            
            ans.mat=matrix(NA_real_, length(vec), R)
            cur.col=0
            for(i in 1:ncol(this.ans)){
                if(remainder[i]==0) next
                tmp=get.combn(vec[-this.ans[,i]], n[-1], remainder[i])
                tmpans=rbind(matrix(rep(vec[this.ans[,i]],ncol(tmp)),n[1]), tmp)
                ans.mat[,cur.col+1:ncol(tmpans)]=tmpans
                cur.col=cur.col+ncol(tmpans)
            }
            ans.mat[,1:cur.col,drop=FALSE]
        }
        ans.mat=get.combn(1:N, n, B)
#    }
    identity=which(colSums(abs(ans.mat-1:N))==0)
    if(length(identity)==0){
        ans.mat[,1]=1:N
    }else{
        ans.mat[,identity]=ans.mat[,1]
        ans.mat[,1]=1:N
    }
    perms=ans.mat[order(trt),]
    attr(perms,'.Random.seed')=save.seed
    attr(perms,'trt')=trt
    perms
}

