dncp=function(obj, ...) UseMethod("dncp")

dncp.parncpt=function(obj, fold=FALSE, ...)
{
    ans=function(x)dnorm(x,obj$mu.ncp,obj$sd.ncp)
    if(fold) fold(ans) else ans
}

dncp.nparncpt=function(obj, fold=FALSE, ...) {
    d.ncp = function(xx) {
        xx = outer(xx, obj$all.mus, "-")
        xx = sweep(xx, 2, obj$all.sigs, "/")
        d = sweep(dnorm(xx), 2, obj$all.sigs, "/")
        drop(d %*% obj$beta)
    }
    if(fold) fold(d.ncp) else d.ncp    
}

dncp.sparncpt=function(obj, fold=FALSE, ...) 
{
    d.ncp=function(x) obj$par*dncp(obj$parfit)(x) + (1-obj$par)*dncp(obj$nparfit)(x)
    if(fold) fold(d.ncp) else d.ncp
}

dncp.nparncpp=function(obj, reflect=TRUE,...)
{
  .NotYetImplemented()
}


dncp.parncpF=function(obj,...)
{
    ans=function(x)dchisq(x/obj$gamma2, obj$data$df1, obj$delta0/obj$gamma2)/obj$gamma2
    ans
}

dncp.nparncpF=function(obj,...)        # p in the paper
{   ## depends on mus, sigs
    d.ncp=function(xx){
        z=function(k, u) dchisq(u/obj$all.gam2s[k], obj$data$df1, obj$all.mus[k]/obj$all.gam2s[k]-obj$data$df1)/obj$all.gam2s[k]
        qx=outer(1:(length(obj$all.mus)-2), xx,  z)
        drop(obj$beta %*% qx)
    }
    d.ncp
}

dncp.sparncpF=function(obj,...) 
{
    d.ncp=function(x) obj$par*dncp(obj$parfit)(x) + (1-obj$par)*dncp(obj$nparfit)(x)
    d.ncp
}
