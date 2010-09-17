loadOrInstall=function(package, dependencies="Depends",...)
{
    stopifnot(all(is.character(package)))
    ddd=list(...)
    singlePKG=function(pkg)
    {
        load.rslt=do.call('require',c(ddd, package=pkg, character.only=TRUE))
        if(isTRUE(load.rslt)) return(TRUE)
        do.call('install.packages',c(ddd, pkgs=pkg,dependencies=dependencies))
        do.call('require', c(ddd, package=pkg, character.only=TRUE, quietly = TRUE))
    }
    rslt=sapply(package, singlePKG)
    if(any(!rslt)) {
        ans=FALSE
        attr(ans,'failed')=package[!rslt]
    }else ans=TRUE
    ans
}
