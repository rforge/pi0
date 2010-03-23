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
        drop(d %*% ob$beta)
    }
    d.ncp
}

#dncp.sparncp=function(obj,...) dncp.nparncp(obj,...)
