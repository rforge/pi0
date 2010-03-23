mdt=function(obj, ...) UseMethod('mdt') ## marginal density of t

mdt.parncp=function(obj, ...)
{
    single.df.ans=function(x, df) obj$pi0*dt(x, df )+(1-obj$pi0)*dtn.mix(x, df,obj$mu.ncp, obj$sd.ncp, ...)
    df.unique=sort(unique(obj$data$df))
    if(length(df.unique)==1) return(function(x)single.df.ans(x,obj$data$df[1]))

    
    ans=function(x){    # discrete mixure of many distinct df's
        dftab=table(obj$data$df)
        prop=dftabl/sum(dftab)
        sums=numeric(length(x))
        for(i in seq(along=dftab) ) sums=sums+prop[i]*single.df.ans(x, dftab[i])
        sums
    }
    ans
}

mdt.nparncp=function(obj, ...)
{
    single.df.ans=function(x,df)
    {   sums=obj$pi0*dt(x, df)
        for(k in 1:length(obj$beta) {
            tmp=dtn.mix(x,df, obj$all.mus[k], obj$all.sigs[k],...)
            if(any(tmp<0) || any(is.na(tmp))) {
                warning("Noncentral density unreliable. I switched to exact density function")
                tmp=dtn.mix(x,df, obj$all.mus[k], obj$all.sigs[k], approximation='none')
            }
            sums=sums+(1-obj$pi0)*obj$beta[k]*tmp
        }
        sums
    }
    df.unique=sort(unique(obj$data$df))
    if(length(df.unique=1)) return(function(x)single.df.ans(x, obj$data$df[1]))

    ans=function(x){    # discrete mixure of many distinct df's
        dftab=table(obj$data$df)
        prop=dftabl/sum(dftab)
        sums=numeric(length(x))
        for(i in seq(along=dftab) ) sums=sums+prop[i]*single.df.ans(x, dftab[i])
        sums
    }
    ans
}

#mdt.sparncp=function(obj, ...){ ## to be implemented
#    
#}
