##############################################################################################
# 
#  CLIENT      BSSI
# 
#  PROJECT     R&D/GENE_LEVEL_2012/KMR
# 
#  PROGRAM		 pchibarsq_LQ_050312.R
# 
#  SOURCE			 
# 
#  DESCRIPTION	Chi-Bar-Square distribution with nonnegatity cone constraints R (Ref: Shapiro, 1988, Int. Stat. Rev., 56, 49; Kudo, 1963, Biometrika, 50, 403)
#                  a.	pchibarsq: Distribution function of chi-bar-square distribution with nonnegativity cone constraint
#                    i.	q: vector of quantiles, as in pchisq. 
#                    ii.	V: a positive-definite matrix, defining the distance measure used when projecting onto the cone. 
#                    iii.	lower.tail: logical, the same as in pchisq. 
#                    iv.	log.p: logical, the same as in pchisq. 
#                  b.	wchibarsq: Computing the mixing proportions for the chi-bar-square distribution.
#                    i.	V: a positive definite matrix as in pchibarsq.
#                  c.	mchibarsq: Computing the moments of the chi-bar-square distribution. 
#                    i.	V: a positive definite matrix as in pchibarsq.
#                    ii.	order: A positive integer vector of the order of moments to be computed.
#
# 
#  LOCATION		//StorageServer01/clients/EL/METHODS_2011/RKHS/PROGRAMS/DEV/
# 
#  INPUT(0)		
#                 
#  OUTPUT(0)		
#
# PROGRAMMER		LONG QU / MAY 03, 2012
#
# MODIFIER      
#
# QC/DATE		    
#
##############################################################################################

# require(mvtnorm)

pchibarsq=function(q, V, lower.tail=TRUE, log.p=FALSE)
{
  n=nrow(V)
  
  wts=wchibarsq(V) 
    
  ans=pchisq(q, 0, lower.tail=FALSE)*wts[1L]+
      pchisq(q, n, lower.tail=FALSE)*wts[n+1L]
  for(i in seq_len(n-1)){
    ans=ans+pchisq(q, i, lower.tail=FALSE)*wts[i+1L]
  }
  ans[q<=0]=1

  ans = if(isTRUE(lower.tail)) 1-ans else ans
  if(isTRUE(log.p)) log(ans) else ans
}

wchibarsq=function(V) ## weights
{## Ref: Kudo A. 1963. Biometrika, 50, pg 403  (page 414)
  seed=get.seed()
  stopifnot(is.matrix(V) && diff(dim(V))==0L && nrow(V)>0L)
  n=nrow(V)
  P=function(idx){## V, n, i
    pmvnorm(rep(0,n-i), rep(Inf, n-i), sigma=solve(V[-idx,-idx,drop=FALSE]))*
    pmvnorm(rep(0, i), rep(Inf, i), 
            sigma=V[idx,idx,drop=FALSE]-V[idx, -idx, drop=FALSE]%*%solve(V[-idx,-idx,drop=FALSE], V[-idx, idx, drop=FALSE])
           )
  }
  ans=numeric(n+1L)
  ans[1]=pmvnorm(rep(0,n),rep(Inf,n),sigma=solve(V))[[1]]
  ans[n+1L]=pmvnorm(rep(0,n), rep(Inf,n), sigma=V)[[1]]
  for(i in safeseq(1L, n-1L, by=1L)) ans[i+1]=sum(combn(n, i, P))
  attr(ans, 'seed')=seed
  ans
}

mchibarsq=function(V, order=1:2) # moments
{
  n=nrow(V)
  wts=wchibarsq(V)
  l2=log(2); lg=lgamma(1:n/2)
  sapply(order,function(ord) weighted.mean(c(0, exp(ord*l2+lgamma(ord+1:n/2)-lg)), wts))
}

if(FALSE){
  safeseq=function(from=1L, to=1L, by=1L,...)
  {
    disc=by*(from-to)
    if(disc>0){
      vector(class(disc), 0L)
    }else seq(from=from, to=to, by=by, ...)
  }
  library(quadprog)
  library(mvtnorm)
  
  V=crossprod(matrix(rnorm(25),5))
  VI=solve(V)
  L=t(chol(V))
  chibarsq=replicate(1e3L, {y=L%*%rnorm(5); -2*solve.QP(VI, VI%*%y, diag(1,5), rep(0,5))$value})
  chibarsq=sort(chibarsq)
  p=pchibarsq(chibarsq, V)
  
  plot(ecdf(chibarsq),new=TRUE)
  lines(chibarsq, p, col=4, lwd=3, lty=3)
  
  mean(chibarsq); mean(chibarsq^2)
  mchibarsq(V)
  
  library(parallel)
  Amat=diag(1,5L); bvec=rep(0,5L)
  cl=makeForkCluster(as.integer(Sys.getenv('OMP_NUM_THREADS')))
  clusterSetRNGStream(cl, 234235)
  chibarsq1e6=unlist(clusterApply(cl, seq_len(1e6L), function(i){y=L%*%rnorm(5); -2*solve.QP(VI, VI%*%y, Amat, bvec)$value})) 
  mean(chibarsq1e6); mean(chibarsq1e6^2)
  mchibarsq(V)
}
