parncpF=function(Fstat, df1, df2, central=TRUE, method=c('L-BFGS-B',"Nelder-Mead", "BFGS", "CG", "SANN"), ...)
{   stopifnot(all(df1>0, df2>0))
     
     method=match.arg(method)
    if       (method=='EM') {
        stop("EM algorithm not implemented")
#        parncpF.em(Fstat,df,central,...)
    }else if (method=='NR') {
        stop("Newton-Raphson algorithm not implemented")
#        parncpF.nr(Fstat,df,central,...)
    }else if (method=='L-BFGS-B') {
        if(central) parncpF.lbfgsb.central(Fstat,df1,df2,...) else parncpF.lbfgsb.noncentral(Fstat,df1,df2,...)
    }else if (method%in%c("Nelder-Mead", "BFGS", "CG", "SANN")) {
        if(central) parncpF.unconstrained.central(Fstat,df1,df2,method=method,,...) else parncpF.unconstrained.noncentral(Fstat,df1,df2,method=method,,...)
    }
}


parncpF.lbfgsb.noncentral=function(Fstat,df1,df2,starts, grids, approximation='none',...)
{
    G=max(c(length(Fstat),length(df1), length(df2)))

    dF.null=df(Fstat,df1,df2)
    obj=function(parms){
            pi0=parms[1]; delta0=parms[2]; gamma2=parms[3]; scale.fact=(1+gamma2)
            Lik=pi0*dF.null+(1-pi0)*dFsnc.mix(Fstat,df1, df2,delta0,gamma2,FALSE,approximation)
#            Lik=pi0*dF.null+(1-pi0)*dt.int(Fstat/scale.fact,df,delta0/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dF.null+(1-pi0)*dFsnc.mix(Fstat,df1, df2,delta0,gamma2,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.noncentral=function(parms) { ### FIXME: analytical derivative not implemented yet
        #loadOrInstall("numDeriv")
        grad(obj, parms)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, 1e-3), upper=c(1-1e-3, 1, 4), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','delta0','gamma2')
    optimFit=optim(starts,obj,gr=deriv.noncentral, method='L-BFGS-B',lower=c(0,0,0),upper=c(1,Inf,Inf),hessian=TRUE,...)
#    optimFit=optim(starts,obj,gr=NULL, method='L-BFGS-B',lower=c(0,0,0),upper=c(1,Inf,Inf),hessian=TRUE,...)

    ll=-optimFit$value
    attr(ll,'df')=3
    class(ll)='logLik'

    ans=list(pi0=optimFit$par[1], mu.ncp=optimFit$par[2]+df1*optimFit$par[3], sd.ncp=sqrt(2*optimFit$par[3]^2*(2*optimFit$par[2]/optimFit$par[3]+df1)), 
            delta0=optimFit$par[2], gamma2=optimFit$par[3], data=list(Fstat=Fstat, df1=df1, df2=df2), 
             logLik=ll, enp=3, par=optimFit$par,
             obj=obj, gradiant=deriv.noncentral(optimFit$par), hessian=optimFit$hessian)
    class(ans)=c('parncpF','ncpest')
    ans
}
parncpF.unconstrained.noncentral=function(Fstat,df1,df2,starts, grids, method='Nelder-Mead', approximation='none',...)
{
    G=max(c(length(Fstat),length(df1), length(df2)))

    logit=make.link('logit')$linkfun
    logit.inv=make.link('logit')$linkinv

    dF.null=df(Fstat,df1,df2)
    obj=function(parms){
            pi0=logit.inv(parms[1]); delta0=exp(parms[2]); gamma2=exp(parms[3]); scale.fact=(1+gamma2)
            Lik=pi0*dF.null+(1-pi0)*dFsnc.mix(Fstat,df1, df2,delta0,gamma2,FALSE,approximation)
#            Lik=pi0*dF.null+(1-pi0)*dt.int(Fstat/scale.fact,df,delta0/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dF.null+(1-pi0)*dFsnc.mix(Fstat,df1, df2,delta0,gamma2,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.noncentral=function(parms) { ### FIXME: analytical derivative not implemented yet
        #loadOrInstall("numDeriv")
        grad(obj, parms)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(-7, -2, -7), upper=c(7, 2, log(4)), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','delta0','gamma2')
    optimFit=optim(starts,obj,gr=deriv.noncentral, method=method,hessian=FALSE,...)
#    optimFit=optim(starts,obj,gr=NULL, method='L-BFGS-B',lower=c(0,0,0),upper=c(1,Inf,Inf),hessian=TRUE,...)

    parncpF.lbfgsb.noncentral(Fstat,df1,df2,starts=c(logit.inv(optimFit$par[1]), exp(optimFit$par[2]), exp(optimFit$par[3])), approximation=approximation,...)
}


parncpF.lbfgsb.central=function(Fstat,df1, df2, starts, grids, approximation='none',...)
{
    G=max(c(length(Fstat),length(df1), length(df2)))

    dF.null=df(Fstat,df1,df2)

    obj=function(parms){
                pi0=parms[1]; delta0=0; gamma2=parms[2]; scale.fact=1+gamma2
                Lik=pi0*dF.null+(1-pi0)*df(Fstat/scale.fact,df1,df2)/scale.fact
                -sum(log(Lik))
    }
    deriv.central=function(parms){### FIXME: analytical derivative not implemented yet
        #loadOrInstall("numDeriv")
        grad(obj, parms)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3), upper=c(1-1e-3, 4), ngrid=c(20,200))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','gamma2')
    optimFit=optim(starts,obj,gr=deriv.central, method='L-BFGS-B', lower=c(0,0),upper=c(1,Inf), hessian=TRUE,...)
#    optimFit=optim(starts,obj,gr=NULL, method='L-BFGS-B', lower=c(0,0),upper=c(1,Inf), hessian=TRUE,...)

    ll=-optimFit$value
    attr(ll,'df')=2
    class(ll)='logLik'
    
#    loadOrInstall("numDeriv")
#    tmp=make.link('logit'); logit=tmp$linkfun; logitinv=tmp$linkinv; dlogitinv=tmp$mu.eta
#    obj.nobound=function(par)obj(c(logitinv(par[1]),exp(par[2])))
#    app.hess.nobound=hessian(obj.nobound, c(logit(optimFit$par[1]), log(optimFit$par[2])))  ## need to consider hitting boundaries
#    app.hess=app.hess.nobound/tcrossprod(c(dlogitinv(logit(optimFit$par[1])), optimFit$par[2]))
#
    ans=list(pi0=optimFit$par[1], mu.ncp=df1*optimFit$par[2], sd.ncp=sqrt(2*optimFit$par[2]^2*df1), 
            delta0=0, gamma2=optimFit$par[2], data=list(Fstat=Fstat, df1=df1, df2=df2), 
             logLik=ll, enp=2, par=optimFit$par,
             obj=obj, gradiant=deriv.central(optimFit$par), hessian=optimFit$hessian)
    class(ans)=c('parncpF','ncpest')
    ans
}

parncpF.unconstrained.central=function(Fstat,df1, df2, starts, grids, method='Nelder-Mead', approximation='none',...)
{
    G=max(c(length(Fstat),length(df1), length(df2)))

    logit=make.link('logit')$linkfun
    logit.inv=make.link('logit')$linkinv
    dF.null=df(Fstat,df1,df2)

    obj=function(parms){
                pi0=logit.inv(parms[1]); delta0=0; gamma2=exp(parms[2]); scale.fact=1+gamma2
                Lik=pi0*dF.null+(1-pi0)*df(Fstat/scale.fact,df1,df2)/scale.fact
                -sum(log(Lik))
    }
    deriv.central=function(parms){### FIXME: analytical derivative not implemented yet
        #loadOrInstall("numDeriv")
        grad(obj, parms)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(-7, -7), upper=c(7, log(4)), ngrid=c(20,200))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','gamma2')
    optimFit=optim(starts,obj,gr=deriv.central, method=method, hessian=FALSE,...)
#    optimFit=optim(starts,obj,gr=NULL, method='L-BFGS-B', lower=c(0,0),upper=c(1,Inf), hessian=TRUE,...)

    parncpF.lbfgsb.central(Fstat,df1,df2,starts=c(logit.inv(optimFit$par[1]), exp(optimFit$par[2])), approximation=approximation,...)
}

## vcov and logLik are defined for the ncpest class in nparncpF.R
#vcov.parncpF=function(obj)
#{
#    obj$hessian
#}
#logLik.parncpF=function(obj)
#{
#    obj$logLik
#}
#coef.parncpF=coefficients.parncpF=function(obj)
#{
#    obj$par
#}
fitted.parncpF=fitted.values.parncpF=function(object, ...)
{
    object$pi0*df(object$data$Fstat, object$data$df1, object$data$df2)+
    (1-object$pi0)*dFsnc.mix(object$data$Fstat, object$data$df1, object$data$df2, object$delta0,object$gamma2,FALSE,...)
}
summary.parncpF=function(object,...)
{
    cat("pi0 (proportion of null hypotheses) =", object$pi0, fill=TRUE)
    cat('mu.ncp=', object$mu.ncp, fill=TRUE)
    cat('sd.ncp=', object$sd.ncp, fill=TRUE)
    cat("delta0  =", object$delta0, fill=TRUE)
    cat("gamma2  =", object$gamma2, fill=TRUE)
    invisible(object)
}
print.parncpF=function(x,...)
{
    summary.parncpF(x,...)
}
plot.parncpF=function(x,...)
{
#    x11(width=8, height=4)
    op=par(mfrow=c(1,2))
    hist(x$data$Fstat, pr=TRUE, br=min(c(max(c(20, length(x$data$Fstat)/100)), 200)), xlab='F',main='F-statistics')
    ord=order(x$data$Fstat)
    lines(x$data$Fstat[ord], fitted.parncpF(x)[ord], col='red', lwd=2)
    d.ncp=function(d) dchisq(d/x$gamma2, x$data$df1, x$delta0/x$gamma2)/x$gamma2
    curve(d.ncp, 0, max(x$data$Fstat), 500, xlab='delta', ylab='density',main='noncentrality parameters')
#    abline(v=c(0, x$delta0), lty=1:2)
    par(op)
    invisible(x)
}

#grid.search=function(obj, lower, upper, ngrid, ...)
#{
#    p=max(c(length(lower),length(upper),length(ngrid)))
#    lower=rep(lower, length=p)
#    upper=rep(upper, length=p)
#    ngrid=rep(ngrid, length=p)
#
#    knot.list=list()
#    for(i in 1:p) knot.list[[i]]=seq(from=lower[i], to=upper[i], length=ngrid[i])
#    names(knot.list)=paste(names(lower),names(upper),names(ngrid),sep='.')
#    
#    grids=do.call(expand.grid, knot.list)
#    ans=apply(grids, 1, obj, ...)
#    if(sum(ans==min(ans))>1) warning('multiple minimums found in grid search')
#    return(unlist(grids[which.min(ans),]))
#}
#
