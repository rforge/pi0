#satterth <-
#function(object, ...) UseMethod('satterth')

satterth.varComp <-
function(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
{
# S3 method for object of class varComp. 
# Lmat:   A matrix with each row being a linear combination of fixed effect parameters of interest. 
# Vbet:   vcov of fixed effect parameter estimates
# svd.VLbet: An optional svd object of vcov of Lmat%*%beta
# X:  Optional X design matrix
# K:  Optional kernal matrix
# V:  Optional vcov of response. 

	if(ncol(model.matrix(object, 'X')) == 0L) return(structure(numeric(0L), individual.df=numeric(0L)))
	
  VVC=vcov(object, what='varComp', drop=TRUE)
  if(missing(Vbet) ) Vbet=vcov(object, what='beta', beta.correction=FALSE)
  if(missing(svd.VLbet)){
    eig=eigen(tcrossprod(Lmat%*%Vbet, Lmat), TRUE)
  }else{
    eig=list(values=svd.VLbet$d, vectors=svd.VLbet$u)
  }
  if(missing(X)) X=model.matrix(object, what='fixed')
  if(missing(K)) K=model.matrix(object, what='K')[object$parms>0]
  if(!is.list(K)) K=list(K)
  if(missing(V)) V=vcov(object, what='Y')
  if(!is.matrix(Lmat)) Lmat=as.matrix(Lmat)
  if(ncol(Lmat)!=ncol(X) && nrow(Lmat)==ncol(X)) Lmat=t(Lmat)
  
  q=sum(eig$values>=sqrt(.Machine$double.eps))
  l=crossprod(eig$vectors[,seq_len(q)], Lmat)
  
  nK=length(K)
  K[[nK+1L]]=diag(1, nrow(X))
  
  VIXUl=solve(V, tcrossprod(X%*%Vbet, l))
  
  ders=matrix(NA_real_, q, nK+1L)
  for(i in seq_len(nK+1L))
    ders[, i]=diag(crossprod(VIXUl, K[[i]]%*%VIXUl))
  
  vs=numeric(q)
  for(m in seq_len(q))
    vs[m]=2*eig$values[m]*eig$values[m]/crossprod(ders[m,], VVC%*%ders[m,])
  
  E=sum(vs/(vs-2)*(vs>2))
  ans=if(E>q) 2*E/(E-q) else 0
  attr(ans, 'individual.df')=vs  
  names(ans)='denominator.df'
  ans
}
