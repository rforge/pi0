parncpt2=function(tstat, df,  ...)
{   stopifnot(all(df>0))
    # if(any(is.infinite(df)) && !all(is.infinite(df)) ) {
        # df[is.infinite(df)]=500
    # }
     method=c('constrOptim')
#     method=match.arg(method)
    if       (method=='EM') {
        stop("EM algorithm not implemented")
#        parncpt.em(tstat,df,zeromean,...)
    }else if (method=='NR') {
        stop("Newton-Raphson algorithm not implemented")
#        parncpt.nr(tstat,df,zeromean,...)
    }else if (method=='constrOptim') {
        parncpt2.constrOptim(tstat,df,...)
    }
}


parncpt2.constrOptim=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; pi1=parms[2]; mu1.ncp=parms[3]; sd1.ncp=parms[4]; mu2.ncp=parms[5]; sd2.ncp=parms[6]; 
			 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
            Lik=pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation='none')+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.obj=function(parms) {#FIXME
        pi0=parms[1]; mu.ncp=parms[2]; sd.ncp=parms[3]; scale.fact=sqrt(1+sd.ncp*sd.ncp); s2=scale.fact*scale.fact
        dt.alt=dtn.mix(tstat, df, mu.ncp, sd.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+(1-pi0)*dt.alt)

        der.pi0=sum( (dt.alt-dt.null) / f )  ## correct

        if(all(is.infinite(df))){
            z.std=(tstat-mu.ncp)/scale.fact
            der.mu=-(1-pi0)*sum( dt.alt*z.std/scale.fact / f)
            der.scale=(1-pi0)*sum( dt.alt*(1-z.std*z.std)/scale.fact / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; t2vs2=t2+df*s2
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
            logC=logK2-(df.half+.5)*log(t2/s2+df)- df.half*mu.ncp*mu.ncp/t2vs2
    #        integral.xv1=.C('intTruncNormVec',n=as.integer(G),r=rep(as.integer(df+1),length=G), mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
    #                                    sd=rep(as.double(1),length=G), lower=numeric(G), upper=rep(Inf,length=G), ans=numeric(G),NAOK=TRUE)$ans
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            
            der.mu=-sum((1-pi0)/f/s2*(tstat*exp(logC)/sqrt(t2vs2)*integral.xv1-mu.ncp*dt.alt))    ## correct

            der.scale=-sum((1-pi0)/f        /s2/scale.fact/t2vs2*(#*dhs.ds)
                dt.alt*(s2*df*(t2-s2)+mu.ncp*mu.ncp*t2vs2)-exp(logC)*mu.ncp*tstat*(t2vs2+df*s2)/sqrt(t2vs2)*integral.xv1)
            )
        }
        der.sd=sd.ncp/scale.fact*der.scale
        c(pi0=der.pi0, mu.ncp=der.mu, sd.ncp=der.sd)
    }
	deriv.obj=function(parms)numDeriv::grad(obj, parms)
	
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, -2, 1e-3), upper=c(1-1e-3, 1-1e-3, 2, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
		obj.restricted=function(parms){
			if(sum(parms[1:2])>1) Inf else
				obj(c(parms[1],parms[2],parms[3],parms[4],parms[3]*-1,parms[4]))
		}
        starts=grid.search(obj.restricted, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
	ui=rbind(diag(1,6)[-c(3,5),], rep(-1:0, c(2,4)), c(0,0,-1,0,1,0))
	ci=rep(0,6); ci[5]=-1
    names(starts)=c('pi0','pi1','mu1.ncp','sd1.ncp','mu2.ncp','sd2.ncp')
    optimFit=try(constrOptim(starts,obj,grad=deriv.obj, ui=ui,ci=ci,hessian=FALSE,...))
    # if(class(optimFit)=='try-error'){
        # optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    # }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }
	optimFit$hessian=numDeriv::hessian(obj, optimFit$par)

    ll=-optimFit$value
    attr(ll,'df')=6
	attr(ll,'nobs')=G
    class(ll)='logLik'

	tau=optimFit$par[2L]/(1-optimFit$par[1L])
	mu=tau*optimFit$par[3L]+(1-tau)*optimFit$par[5L]
	Var=tau*((optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)+
	(1-tau)*((optimFit$par[5L]-mu)^2+optimFit$par[6L]^2)
	
    ans=list(pi0=optimFit$par[1], mu.ncp=mu, sd.ncp=sqrt(Var), tau.ncp=tau, pi1=optimFit$par[2], mu1.ncp=optimFit$par[3], sd1.ncp=optimFit$par[4], 
		mu2.ncp=optimFit$par[5], sd2.ncp=optimFit$par[6], 
		data=list(tstat=tstat, df=df), 
             logLik=ll, enp=6, par=optimFit$par,
             obj=obj, gradiant=numDeriv::grad(obj, optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt2','parncpt','ncpest')
    ans
}

fitted.parncpt2=#fitted.values.parncpt=
function(object, ...)
{
    object$pi0*dt(object$data$tstat, object$data$df)+
    object$pi1*dtn.mix(object$data$tstat, object$data$df,object$mu1.ncp,object$sd1.ncp,FALSE,...)+
	(1-object$pi0-object$pi1)*dtn.mix(object$data$tstat, object$data$df,object$mu2.ncp,object$sd2.ncp,FALSE,...)
}

summary.parncpt2=function(object,...)
{
    cat("pi0 (proportion of null hypotheses) =", object$pi0, fill=TRUE)
    cat("mu.ncp (mean of noncentrality parameters) =", object$mu.ncp, fill=TRUE)
    cat("sd.ncp (SD of noncentrality parameters) =", object$sd.ncp, fill=TRUE)
	 cat(sprintf("tau.ncp=%.5f; mu1.ncp=%.3f; sd1.ncp=%.2f; mu2.ncp=%.3f; sd2.ncp=%.2f",object$tau,object$mu1.ncp,object$sd1.ncp,object$mu2.ncp,object$sd2.ncp), filel=TRUE)
    invisible(object)
}
print.parncpt2=function(x,...)
{
    summary.parncpt2(x,...)
}
plot.parncpt2=function(x,...)
{
	tau=x$tau.ncp
	mu=x$mu.ncp
#    dev.new(width=8, height=4)
    op=par(mfrow=c(1,2))
    hist(x$data$tstat, pr=TRUE, br=min(c(max(c(20, length(x$data$tstat)/100)), 200)), xlab='t',main='t-statistics')
    ord=order(x$data$tstat)
    lines(x$data$tstat[ord], fitted.parncpt2(x)[ord], col='red', lwd=2)
    d.ncp=function(d) tau*dnorm(d, x$mu1.ncp, x$sd1.ncp)+(1-tau)*dnorm(d, x$mu2.ncp, x$sd2.ncp)
    curve(d.ncp, min(x$data$tstat), max(x$data$tstat), 500, xlab=expression(delta), ylab='density',main='noncentrality parameters')
    abline(v=c(0, mu), lty=1:2)
    par(op)
    invisible(x)
}
