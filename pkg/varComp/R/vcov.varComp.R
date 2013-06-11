vcov.varComp <-
function(object, what=c('fixed','beta','random','varComp','var.ratio','tau','response','Y'), drop=TRUE, ...)
{
# S3 method for returning variance-covariance matrix of interesting quantities. 
# what: 'fixed' or 'beta': vcov for fixed effect parameter estimates
#       'random' or 'varComp': vcov for REML estimates
#       'var.ratio' or 'tau': vcov for ratios of variances
#       'response' or 'Y': vcov of response variable. 
  what=match.arg(what)
  what=switch(what, fixed='beta', random='varComp', var.ratio='tau', response='Y', what)
  
  if(what=='Y'){
    K=model.matrix(object, what='K')
    V=Reduce('+', mapply('*', object$varComps, K, SIMPLIFY=FALSE))+diag(object$sigma2, nrow(K[[1]]))
    return(V)    
  }else if(what=='varComp'){
    if(isTRUE(drop)) idx=which(object$parms>0) else idx=seq_along(object$parms)
  
    this.Q2Y=object$residual.contrast
    k=object$working.cor[idx]
    nk=length(k)
    k[[nk+1L]]=diag(1,length(this.Q2Y))
    V=Reduce('+', mapply('*', c(object$parms[idx],1)*object$sigma2, k, SIMPLIFY=FALSE))
    # LI=solve(t(chol(V)))
	LI = t(backsolve(chol(V), diag(1, nrow(V))))
    lik=lapply(k, function(kk) tcrossprod(LI%*%kk, LI))
    liy=LI%*%this.Q2Y
    
    ans=matrix(NA_real_, nk+1L, nk+1L)
    for(i in seq_len(nk+1)){
      for(j in i:(nk+1)){
        tmp=lik[[i]]%*%lik[[j]]
        ans[i,j]=ans[j,i]=-crossprod(liy, tmp%*%liy)+.5*sum(diag(tmp))
      }
    }
    ans=-solve(ans)
    #attr(ans, 'V')=V
    ans
  }else if(what=='tau'){
    if(isTRUE(drop)) idx=which(object$parms>0) else idx=seq_along(object$parms)
    ans=-solve(object$hessian[idx,idx, drop=FALSE])
    ans 
  }else if(what=='beta'){
    X=model.matrix(object, what='X')
    ginv(crossprod(X, solve(vcov(object, what='Y'),X)))
  }
}
