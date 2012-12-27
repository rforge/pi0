FR2dec=function(FR){
    drop(gmp::crossprod(factorialZ(seq_along(FR)-1),rev(FR)))
}
if(FALSE){
dec2FR=function(dec,N){	## old slow but correct implementation
    if(missing(N)){
        k=1L
        fact=1L
        repeat{
            fact=fact*k
            if(fact>dec)break
            k=k+1L
        }
        N=k
    }
    ans=integer(N)
    for(i in 1:N){############ algorithm in Chapter 10 of .NET Test Automation Recipes
        ans[N-i+1]=as.integer(mod.bigz(dec, i))  
        dec=divq.bigz(dec, i)
    }
    ans
}
}

dec2FR=function(dec,N){		## new implementation
    if(missing(N)){
        k=1L
        fact=as.bigz(1L)
        repeat{
            fact=fact*k
            if(fact>dec)break
            k=k+1L
        }
        N=k
    }
	dec
	as.integer(mod.bigz(divq.bigz(dec, factorialZ(N:1-1L)), N:1))
}

if(FALSE){
FR2permvec=function(FR,base=1L){  ## TODO: optimize this function
    stopifnot(base%in%(0:1))
    ans=integer(length(FR))
    cand=seq_along(FR)-(1L-base)
    for(i in 1:length(FR)){ ## see http://en.wikipedia.org/wiki/Factoradic
        idx=FR[i]+1L
        ans[i]=cand[idx]
        cand=cand[-idx]
    }
    ans
}
}

FR2permvec=function(FR, base=1L){  ## C implementation
	if(!is.integer(FR)) FR=as.integer(FR)
	if(!is.integer(base)) base=as.integer(base)
	.Call('FR2permvec', FR, base)
}

######## this is a naive inversion algorithm to check correctness
#permvec2FR0=function(permvec){
#    permvec=permvec-min(permvec)
#    N=length(permvec)
#    ans=cand=seq_len(N)-1
#
#    for(i in 1:N){
#        tmp=which(cand==permvec[i])
#        ans[i]=tmp-1
#        cand=cand[-tmp]
#    }
#    ans
#}
#permvec2FR1=function(permvec,base){## an algorithm that is (often) faster
#    if(missing(base))permvec=permvec-min(permvec)
#    permvec-sapply(seq_along(permvec),function(k)sum(permvec[1:k]<permvec[k]))
#}
permvec2FR=function(permvec){## an equivalent algorithm at http://www.mathe2.uni-bayreuth.de/frib/KERBER/h00/node30.html
#    if(missing(base))permvec=permvec-min(permvec)
    N=length(permvec)
    sapply(seq_len(N),function(k)sum(permvec[k:N]<permvec[k]))
}


dec2permvec=function(dec,N,base=1L) FR2permvec(dec2FR(dec,N),base)
#permvec2dec0=function(permvec) FR2dec(permvec2FR0(permvec))
permvec2dec=function(permvec) FR2dec(permvec2FR(permvec))
