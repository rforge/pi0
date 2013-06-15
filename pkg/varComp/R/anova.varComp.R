anova.varComp <-
function(object, Lmat, alpha=.05, test='Satterthwaite', ...)
{######### FIXME: add KR and LRT
## S3 method for reporting fixed effect parameters from class varComp. 
# object: varComp object
# Lmat: The same as in satterth, default to all parameters (i.e., identity Lmat)
# alpha: Type I error level requested. 
	if(is.null(test[1L])) return(coef(object, 'fixed'))
	test=match.arg(test)

	if('Y'%in%names(object)){
		Y = object$Y
	}else if ('model'%in%names(object)){
		Y = model.response(object$model)
	}else stop("response variable is not recored.")

  
  X=model.matrix(object, what='fixed')
  K=model.matrix(object, what='K')[object$parms>0]
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat)) Lmat=as.matrix(Lmat)
  if(ncol(Lmat)!=ncol(X) && nrow(Lmat)==ncol(X)) Lmat=t(Lmat)
  
  this.V=vcov(object, what='Y')
  this.Vbet=solve(crossprod(X, solve(this.V,X)))
  this.bet=this.Vbet%*%crossprod(X, solve(this.V, Y))
  
  est=drop(Lmat%*%this.bet)
  LVL=tcrossprod(Lmat%*%this.Vbet, Lmat)
  svd.LVL=svd(LVL)
  LVLI=tcrossprod(sweep(svd.LVL$v, 2, ifelse(svd.LVL$d>sqrt(.Machine$double.eps), 1/svd.LVL$d, 0), '*'), svd.LVL$u) 
  rk=sum(svd.LVL$d>sqrt(.Machine$double.eps))
  
  ses=sqrt(diag(LVL))
  t.dfs=numeric(nrow(Lmat))
  for(i in seq_len(nrow(Lmat))) 
    t.dfs[i]=satterth(object=object, Lmat[i,,drop=FALSE], Vbet=this.Vbet, 
                      #svd.VLbet=list(d=svd.LVL$d[i], u=svd.LVL$u[i,,drop=FALSE], v=svd.LVL$v[i,,drop=FALSE]),
                      X=X, K=K, V=this.V )
#  t.dfs=apply(Lmat, 1L, satterth, object=object, Vbet=this.Vbet, X, K=K, V=this.V)  ## X has naming conflict
  ll=est - ses*qt(1-alpha/2, t.dfs)
  ul=est + ses*qt(1-alpha/2, t.dfs)
  tstats=est/ses
  t.pval=pt(abs(tstats), t.dfs, lower.tail=FALSE)*2
  #debug(satterth.varComp)
  Fstat=crossprod(est, LVLI%*%est)/rk
  F.df=satterth(object, Lmat=Lmat, Vbet=this.Vbet, svd.VLbet=svd.LVL, K=K, V=this.V, X)
  F.p=pf(Fstat, rk, F.df, lower.tail=FALSE)
  
  ans=cbind(estimate=est, se=ses, lower=ll, upper=ul, t=tstats, df=t.dfs, `Pr(>|t|)`=t.pval)
  attr(ans, 'Overall')=cbind(F=drop(Fstat), df1=rk, df2=F.df, `Pr(>F)`=drop(F.p))
  rownames(ans)=rownames(Lmat)
  ans
}
