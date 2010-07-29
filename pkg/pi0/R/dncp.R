dncp=function(obj, ...) UseMethod("dncp")

dncp.parncp=function(obj, ...)
{
    ans=function(x)dnorm(x,obj$mu.ncp,obj$sd.ncp)
    ans
}

dncp.nparncp=function(obj, ...) {
    d.ncp = function(xx) {
        xx = outer(xx, obj$all.mus, "-")
        xx = sweep(xx, 2, obj$all.sigs, "/")
        d = sweep(dnorm(xx), 2, obj$all.sigs, "/")
        drop(d %*% obj$beta)
    }
    d.ncp
}

dncp.sparncp=function(obj,...) 
{
    d.ncp=function(x) obj$par*dncp(obj$parfit)(x) + (1-obj$par)*dncp(obj$nparfit)(x)
    d.ncp
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
