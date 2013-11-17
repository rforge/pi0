if (FALSE) {  ## old code
get.perm2.mat <-
function(trt, B=100L) ## trt needs have exactly 2 levels; other situations are not implemented yet
{
    if(!exists('.Random.seed', envir=globalenv())) runif(1L)
    save.seed=get('.Random.seed', envir=globalenv())
    n=table(trt)
    if(length(n) != 2L) stop('"get.perm2.mat" only works when "trt" has exactly two levels.')
    N=sum(n); m=min(n)
    #perms=replicate(B, sort(sample(N, m)))
    perms=combn2R(N, m, R=B, sample.method="noReplace")
    idx1=sort(which(as.factor(trt)==names(n[which.min(n)])))
    colid=which(apply(perms == idx1 , 2L, all))
    if (length(colid)==0L) {
        perms[,1L]=idx1
    }else if(colid!=1L){
        perms[, colid] = perms[,1L]
        perms[,1L] = idx1 
    }
    attr(perms,'.Random.seed')=save.seed
    perms
}


permuteTrt <-
function(trt, B=100L) ## permutation matrix for one way design
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
            if(FALSE){  ## this is old implementation that is buggy, but still useful
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
                    this.ans=combn2R(length(vec), n[1],R=R) ## BUG: this does not allow duplicates
                    remainder=rep(1L, R)
                }
            }else{
                remainder=table(bdSample(seq(choose(length(vec), n[1L])), 
                                         R, 
                                         get.choose(length(vec)-n[1L], n[-1L])
                                         )
                                )
                this.ans=combn(length(vec), n[1L])[,as.integer(names(remainder)),drop=FALSE]
            }
            
            ans.mat=matrix(NA_real_, length(vec), R)
            cur.col=0
            for(i in 1:ncol(this.ans)){
                if(remainder[i]==0) next
                tmp=Recall(vec[-this.ans[,i]], n[-1], remainder[i])
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

bdSample=function(x, n, ubounds)
# sample n elements of vector x with replacement subject to the number of each element being sampled to be no larger than the corresponding upper bound ubounds
{
    N=length(x)
    x0=seq(N)
    if(length(ubounds)!=N) ubounds=rep(ubounds, length=N)
    stopifnot(sum(ubounds)>n)
    counts=integer(N)
    idx=rep(TRUE, N)
    m=n
    ans=c()
    while(m>0L){
        tmp=x0[idx][sample(sum(idx), m, replace=TRUE)]
        tmpCt=table(tmp); unq=as.integer(names(tmpCt))
        outOfBd= counts[unq]+tmpCt > ubounds[unq]
        if(any(outOfBd)){
            tmpCt[outOfBd]=ubounds[unq[outOfBd]] - counts[unq[outOfBd]]
            idx[unq[outOfBd]]=FALSE
        }
        ans=c(ans, rep(unq, tmpCt))
        counts[unq]=counts[unq]+tmpCt
        m=m-sum(tmpCt)
    }
    if(length(ans)==1) x[ans] else x[sample(ans)]
}

mchooseZ=function(N, n)  ## multinomial coef for a vector of n
{
   ans = chooseZ(N, n[1L])
   if(length(n)==1L) ans else ans*Recall(N-n[1L], n[-1L])
}
}

nparts=function(n)
  factorialZ(sum(n))/prod(c(factorialZ(n), factorialZ(table(n))))  ## total number of distinct trt assignments 

permuteTrt <-
function(trt, B=100L, idxOnly = FALSE) ## permutation matrices for one way design
## returns a length(trt) by B matrix, if B is not larger than multinomial coefficient. 
{
    n=table(trt)
    cn=cumsum(n)
    ntrts=length(n)
    N=cn[ntrts]
    #mc=mchooseZ(N, n)
    ordn=order(n, decreasing=TRUE)

    SP=nparts(n)  
    part0=split(seq_len(N),trt); 
    
    if(B>=SP){ # list all partitions
        sp=setparts(n)
        B=ncol(sp)
#        for(i in seq(ntrts))  ans[[i]]=matrix( apply(sp==i,2L,function(xx)sort(which(xx)) ), ncol=B)
        ans=split(row(sp),sp); ## this line and the next replace the previous line
        for(i in seq(ntrts)) dim(ans[[i]])=c(length(ans[[i]])/B,B)   
        names(ans)=names(n)[ordn]
        
        # swapping the original assignment to the first permutation        
        flag=TRUE
        for(b in seq(B))  
            if(setequal(part0, lapply(ans, '[', , b))){ idx=b; flag=FALSE; break }
        if(isTRUE(flag)){
          warning("The first permutaiton may not be the original assignment.")
        }
        for(i in seq(ntrts)) {tmp=ans[[i]][,b]; ans[[i]][,b]=ans[[i]][,1L]; ans[[i]][,1L]=tmp} 

        if(isTRUE(idxOnly)){
			warning("'idxOnly=TRUE' has not been implemented yet. Full results are returned.")			
        }#else{
            attr(ans, 'idx') = NA_character_
    		class(ans)='permutedTrt'
        #}
    }else{   #sample from all permutations using factoradic number. Ideally, a sample from 1:SP should work, but how to do this without enumerating all SP possibilities using setparts?

        decfr=HSEL.bigz(factorialZ(N), B)
        idx=which(decfr==0L)
        if(length(idx)>0L) decfr[idx]=decfr[1L]
        decfr[1L]=as.bigz(0L)
		
		decfrCC=drop(as.character(decfr))  ## for speed only
		if(isTRUE(idxOnly)) {	## save memory
			ans = lapply(part0, as.matrix)
			attr(ans, 'idx') = decfrCC
			class(ans) = 'permutedTrt'
			return(ans)
		}
		
       ans=lapply(sapply(part0,length), matrix, data=NA_integer_, ncol=B)
		buff = integer(N); buff[1L]
       for(b in seq(B)){
            # perm=dec2permvec(decfr[b],N)  ## This subsetting decfr[b] is the slowest part!
            perm=dec2permvec(decfrCC[b],N)  ## This change speeds up for about 8~9X.
            # for(i in seq(ntrts)) ans[[i]][,b]=sort.int(perm[part0[[i]]])
			 for(i in seq(ntrts)) ans[[i]][,b]=.Call(radixSort_prealloc, perm[part0[[i]]], buff)  ## radix sort with pre-allocated buffer space
       }
       names(ans)=names(part0)
		attr(ans, 'idx') = NA_character_
		class(ans)='permutedTrt'
    }
    ans
}

nperms.permutedTrt=function(permutedTrt)
{
	if(is.na(attr(permutedTrt, 'idx')[1L])) {
			ncol(permutedTrt[[1L]])
	}else   length(attr(permutedTrt, 'idx'))
}

ntrt.permutedTrt=function(permutedTrt)
{
	sapply(permutedTrt, nrow)
}

trt.permutedTrt=function(permutedTrt)
{
	ans=factor(rep(1L, sum(sapply(permutedTrt,nrow))), levels=names(permutedTrt))
	for(i in seq_along(permutedTrt)) ans[permutedTrt[[i]][,1L]] = i
	ans
}
