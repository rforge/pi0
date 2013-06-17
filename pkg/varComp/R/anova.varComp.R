anova.varComp <-
function (object, ..., test = c("KR", 'Satterthwaite'), L)
{
	cl=match.call(expand.dots=FALSE)
	if('...'%in%names(cl)){
		ddd.names=sapply(cl[['...']], deparse)
		ddd = list(...)
		nObj=length(ddd)
		nObj=nObj + 1L
		ddd[[nObj]] = object
		names(ddd) = c(ddd.names, deparse(substitute(object)))
		if(!all(sapply(ddd, inherits, what='varComp'))) stop("All `...` objects need to inherit the `varComp` class.")
		if(diff(range(sapply(ddd, nobs))) != 0 ) stop("The number of observations needs to be the same in all model fits.")
		n=nobs(ddd[[1L]])
		
		Xmats = lapply(ddd,  model.matrix, what='fixed')
		qrs = lapply(Xmats, qr)
		rks = sapply(qrs, '[[', 'rank')
		ork=order(rks)
		rks=rks[ork]
		
		ddd=ddd[ork]
		Xmats=Xmas[ork]
		qrs=qrs[ork]
		# (`raw.F`=drop(Fstat), `scale.F` = scale.overall, df1=rk, df2=F.ddf, `Pr(>F)`=drop(F.p))
		
		nas=rep(NA_real_, nObjs)
		ans = data.frame(`raw.F`=nas, `scale.F` = nas, df1=nas, df2=nas, `Pr(>F)` = nas, check.names=FALSE)
		
		if('model'%in%names(ddd[[1L]])){
			hasInt = (attr(ddd[[1L]]$model, 'intercept') == 1)
		}else{
			ones = rep(1, n)
			hasInt = ( mean((tcrossprod(qr.Q(qrs[[1]]))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}
		diffLmat=function(qr0, X1)
		{
			z=(diag(1, nrow(X1))- tcrossprod(qr.Q(qr0)))%*%X1
			qrz=qr(z)
			 pivot <- attr(qrz, "pivot") and oo <- order(pivot)
			 qr.R(qrz)[seq_len(qrz$rank), oo, drop=FALSE]
		}
		if(hasInt){
			Lmat = diffLmat(qr(rep(1,n)), Xmats[[1L]])
		}else{
			Lmat = diag(1, ncol(Xmats[[1L]]))
		}
		for(m in seq_len(nObj)){
			tmp = fixef(ddd[[m]], Lmat = Lmat, test=test)
			ans[m, ] = attr(attr(tmp, 'anova'), 'Overall')
			if(m<=nObj)	Lmat = diffLmat(qrs[[m]], Xmats[[m+1L]])
		}
		return(ans)
	}
	
	### one ojbect is given
	if(!missing(L)){
		tmp = fixef(object, Lmat=L, test=test)
	}else tmp = fixef(object, test=test)
	tmp.aov=attr(tmp, 'anova'
	ans = data.frame(`raw.F` = tmp.aova[,'raw.t']^2,
					 `scale.F` = tmp.anova[,'scale.t']^2, 
					 `df1` = 1, 
					 `df2` = tmp.anova[,'df'], 
					 `Pr(>F)` = tmp.anova[, 'Pr(>|t|)']
					)
	if(missing(L)){
		warning("Current implementation will test each fixed effect parameter separately when only one `varComp` object is provided.")
		X=model.matrix(object)
		n = nobs(X)
		if('model'%in%names(object)){
			hasInt = (attr(object$model, 'intercept') == 1)
		}else{
			ones = rep(1, n)		
			hasInt = ( mean((tcrossprod(qr.Q(qr(X)))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}	
		if(hasInt){
			tmp = fixef(object, Lmat = diffLmat(qr(rep(1,n)), X), test=test)
			tmp.aov=attr(tmp, 'anova')
		}
	}
	ans = rbind(ans, attr(tmp.aov, 'Overall'))
	rownames(ans) = c(rownames(L), 'Overall')
	return(ans)

}
