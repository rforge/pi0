loadOrInstall=function(package, dependencies="Depends",...)
{
    for.one.package=function(pkg, ...)
    {
        load.rslt=suppressWarnings(require(pkg, character.only=TRUE,...))
        if(isTRUE(load.rslt)) return(TRUE)
        install.packages(pkg,dependencies=dependencies,...)
        require(pkg, character.only=TRUE, quietly = TRUE, ...)
    }
    rslt=sapply(package, for.one.package)
    if(any(!rslt)) {
        ans=FALSE
        attr(ans,'failed')=package[!rslt]
    }else ans=TRUE
    ans
}
