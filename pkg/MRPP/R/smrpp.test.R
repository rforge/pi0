if (FALSE) {
### This is obsoleted code that relies on Witten & Tibshirani (JASA, 2010) criterion. 

smrpp.test.default <-
function(y, trt, B=nparts(table(trt)), permutedTrt, wtmethod=0, eps=1e-8, spar=seq(1,sqrt(ncol(y)),length=100L), ...) ## this uses C code
## y is a dist object; wtmethod: 0=sample size-1; 1=sample size
{
    p1=2  ## order of univariate minkowski dist, to the power of p1
    p2=2  ## order of multivairate minkowski dist
    if(missing(y) ) stop('dist object missing or incorrect')
    if(!is.matrix(y) && !is.data.frame(y)) y= as.matrix(y)
    R=ncol(y)
    N= nrow(y)
    spar=sort(spar[spar>=1 & spar<=sqrt(R)])
    wtmethod=as.numeric(wtmethod[1])
    if(missing(trt)) {  ## recoving trt from the first permutation
      trt=integer(N)
      for(i in seq_along(permutedTrt)) trt[permutedTrt[[i]][,1L]]=i
    }
    if(missing(permutedTrt)) {
        permutedTrt=permuteTrt(trt,B)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
    B=ncol(permutedTrt[[1]])
    #if(missing(cperm.mat)) cperm.mat=apply(permutedTrt, 2, function(kk)(1:N)[-kk])
    tabtrt=table(trt)
    ntrt=length(tabtrt)
    
    get.bcss=function(thisPerm){
      ans=numeric(R)
      for(r in seq(R)){
        thisDist=dist(y[,r,drop=FALSE], method='minkowski', p=p1)^p1
        ans[r]=sum(thisDist)*2/(N-wtmethod) - .Call('mrppstats', thisDist, thisPerm, wtmethod, PACKAGE='MRPP')
      }
      ans
    }
    sthresh0=function(x,cutoff){  ## this is the general soft-thresholding function
      sign(x)*pmax(abs(x)-cutoff,0)
    }
    sthresh=function(x,cutoff){  ## this is the soft-thresholding function used here when x and cutoff are both non-negative
      pmax(x-cutoff,0)
    }
    get.wt=function(bcss, s){
      ap=pmax(bcss,0)
      tmp=sthresh(ap,0)
      w=tmp/sqrt(sum(tmp*tmp))
      if (sum(w) <= s) return(w)
      s2=s*s
      obj=function(delta){
        tmp=sthresh(ap, delta)
        sum(tmp)^2 - s2*sum(tmp*tmp)
      }
      uap=unique(ap); tol=sqrt(.Machine$double.eps)
      d=uniroot(obj, c(0, max(uap[-which.max(uap)], max(uap)-tol)))$root
      tmp=sthresh(ap,d)
      tmp/sqrt(sum(tmp*tmp))
    }
    ycol=col(y)
    if(p2==2){
      get.wdist=function(wt,p2){
        #wy=sweep(y,2L, sqrt(wt), '*', check.margin=FALSE)
        wy=y*sqrt(wt)[ycol]   ## this is much faster
        dist(wy)
      }
    }else {
      get.wdist=function(wt,p2){
        #wy=sweep(y,2L,sqrt(wt),'*')
        wy=y*sqrt(wt)[ycol]   ## this is much faster
        dist(wy, method='minkowski', p=p2)
      }
    }    
    stats=numeric(B)
    for(b in seq(B)){
      thisPerm=lapply(permutedTrt, '[', , b, drop=FALSE)
      bcss=get.bcss(thisPerm)
      
      wmrpp.p=numeric(length(spar))
      for(s.i in seq_along(spar)){
        wt=get.wt(bcss, spar[s.i])
        wdist=get.wdist(wt,p2)
        tmp=.Call('mrppstats',wdist, permutedTrt, wtmethod, PACKAGE='MRPP')
        wmrpp.p[s.i]=mean(tmp[b]-tmp >=-eps)
      }
      plot(wmrpp.p~spar, main=b)
      stats[b]=min(wmrpp.p)
      if(b==1L){
        s0=spar[which.min(wmrpp.p)]
        wt0=get.wt(bcss, s0)
      }
    }
    
    #stats=.Call('mrppstats',y,permutedTrt, as.numeric(wtmethod[1L]), PACKAGE='MRPP')
    ans=list(statistic=c("Sparse MRPP statistic"=stats[1L]), all.statistics=stats, weights=wt0,
             p.value=mean(stats[1]-stats>=-eps), parameter=c("number of permutations"=B, 'group weight method'=wtmethod[1L], 'Smoothing'=s0, '#selected variables'=sum(wt0>0)),
             data.name=dname, #  .Random.seed=attr(permutedTrt,'.Random.seed'),  ## this is commented out since the random number seed in gmp:::urand.bigz will be used. 
             method=sprintf('%d-sample Sparse MRPP test',ntrt)
             )
    class(ans)='htest'
    ans
}

}

smrpp.defaultSpar <- 
function(dp.dw, max.ratio=2, nspar=100L)
{
# Let k be the ratio of maximum weight to the minimum weight.
# When spar increases, this ratio decreases. 
# Whenever k drops below max.ratio, the corresponding spar will be treated as Inf (with k=1)

    stopifnot(max.ratio>1 && nspar>=3L)
    k=max.ratio; 
    if(!is.matrix(dp.dw) && is.numeric(dp.dw)) {dropf=drop; dp.dw=matrix(dp.dw, 1L)} else dropf=I
    B=nrow(dp.dw); R=ncol(dp.dw)
    ans=matrix(NA_real_, B, nspar)
    for(b in seq(B)){
        m1=min(dp.dw[b,])
        m2=min(dp.dw[b, dp.dw[b,]>m1])
        mR=max(dp.dw[b,])
        min. = (m2-m1)*.5/R*sum(dp.dw[b,]==m1)
        max. = .5*( (k*mR-m1)/(k-1)-mean(dp.dw[b,]))
        ans[b,]=c(10^seq(from=log10(min.), to=log10(max.), length=nspar-1L), Inf)
    }
    unclass(dropf(ans))
}

smrpp.penWt <-
function(dp.dw, spar, simplify=TRUE)
# dp.dw is a vector or matrix of derivative of MRPP p-value to the weights. 
# If dp.dw is a matrix, each row is treated as a permutation. 
# spar is a scalar tuning parameter that controls the influence of the heterogeneity of weights
{
    if(!is.matrix(dp.dw) && is.numeric(dp.dw)) dp.dw=matrix(dp.dw, 1L)
    B=nrow(dp.dw)
    R=ncol(dp.dw)+0.0
    if(is.matrix(spar)) stopifnot(nrow(spar)==nrow(dp.dw))
    else spar=matrix(spar, B, length(spar), byrow=TRUE)
    L=ncol(spar)
    sparMin=smrpp.defaultSpar(dp.dw, nspar=3L)[,1L]
    spar=structure(pmax(spar, sparMin), dim=dim(spar))

    fact=.5/spar
    mean.dp.dw=rowMeans(dp.dw)

#   keep the line below for checking correctness of the C implementation
#    f=function(del) fact[b,l]*sum(pmax(-del-dp.dw[b, ], 0)) - R  ## depends on l and b in the enclosing env
    f=function(del) .Call('objSolveDelta', del, dp.dw, fact, b, l, R, PACKAGE='MRPP')  ## depends on l and b in the enclosing env
    ans=array(NA_real_, dim=c(L, B, R))
    for(b in seq_len(B)){
        rg=-rev(range(dp.dw[b,]))
        for(l in seq_len(L)){
            if(is.infinite(spar[b,l])) {ans[l, b, ]=1; next}
            d=-mean.dp.dw[b]-2*spar[b,l]
            if(d >= rg[1]) d=uniroot(f, rg, tol=1e-15)$root
            ans[l, b, ]=fact[b,l]*pmax(-d-dp.dw[b,], 0)
        }
    }
    if(isTRUE(simplify)) drop(ans) else ans
}

smrpp.test.default <-
function(y, trt, B=nparts(table(trt)), permutedTrt, wtmethod=0, eps=1e-8, spar, ...) ## this uses C code
## y is a dist object; wtmethod: 0=sample size-1; 1=sample size
{
    if(missing(y) ) stop('dist object missing or incorrect')
    if(!is.matrix(y) && !is.data.frame(y)) y= as.matrix(y)
    R=ncol(y)
    N= nrow(y)
    wtmethod=as.numeric(wtmethod[1])
    if(missing(trt)) {  ## recoving trt from the first permutation
      trt=integer(N)
      for(i in seq_along(permutedTrt)) trt[permutedTrt[[i]][,1L]]=i
    }
    if(missing(permutedTrt)) {
        permutedTrt=permuteTrt(trt,B)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
    B=ncol(permutedTrt[[1]])
    #if(missing(cperm.mat)) cperm.mat=apply(permutedTrt, 2, function(kk)(1:N)[-kk])
    tabtrt=table(trt)
    ntrt=length(tabtrt)

    dp.dw=get.dp.dw.kde(y, permutedTrt, test = TRUE, wtmethod=wtmethod)  # B x R matrix
    sparMinMax=apply(dp.dw, 1L, smrpp.defaultSpar, nspar=3L, ...)
    sparMin=min(sparMinMax[1L,])
    if(missing(spar)){
        nspar=100L
        sparMax=max(sparMinMax[2L,])
        spar=c(10^seq(from=log10(sparMin), to=log10(sparMax), length=nspar-1L), Inf)
    } else {
        sparMax=max(spar, sparMinMax[2L,])
        spar=sort(spar[ spar>=sparMin ])
        if(all(is.finite(spar))) spar=c(spar, Inf)
        nspar=length(spar)
    }

#    get.bcss=function(thisPerm){
#      ans=numeric(R)
#      for(r in seq(R)){
#        thisDist=dist(y[,r,drop=FALSE], method='minkowski', p=p1)^p1
#        ans[r]=sum(thisDist)*2/(N-wtmethod) - .Call('mrppstats', thisDist, thisPerm, wtmethod, PACKAGE='MRPP')
#      }
#      ans
#    }
#    sthresh0=function(x,cutoff){  ## this is the general soft-thresholding function
#      sign(x)*pmax(abs(x)-cutoff,0)
#    }
#    sthresh=function(x,cutoff){  ## this is the soft-thresholding function used here when x and cutoff are both non-negative
#      pmax(x-cutoff,0)
#    }
#    get.wt=function(bcss, s){
#      ap=pmax(bcss,0)
#      tmp=sthresh(ap,0)
#      w=tmp/sqrt(sum(tmp*tmp))
#      if (sum(w) <= s) return(w)
#      s2=s*s
#      obj=function(delta){
#        tmp=sthresh(ap, delta)
#        sum(tmp)^2 - s2*sum(tmp*tmp)
#      }
#      uap=unique(ap); tol=sqrt(.Machine$double.eps)
#      d=uniroot(obj, c(0, max(uap[-which.max(uap)], max(uap)-tol)))$root
#      tmp=sthresh(ap,d)
#      tmp/sqrt(sum(tmp*tmp))
#    }

    ycol=col(y)
#    if(p2==2){
      get.wdist=function(wt,p2){ # p2 is a placeholder only for potential extensions
        #wy=sweep(y,2L, sqrt(wt), '*', check.margin=FALSE)
        wy=y*sqrt(wt)[ycol]   ## this is much faster
        dist(wy)
      }
#    }else {
#      get.wdist=function(wt,p2){
#        #wy=sweep(y,2L,sqrt(wt),'*')
#        wy=y*sqrt(wt)[ycol]   ## this is much faster
#        dist(wy, method='minkowski', p=p2)
#      }
#    }
    
    stats=numeric(B)
    all.weights=smrpp.penWt(dp.dw, spar, simplify=FALSE)

    for(b in seq(B)){
#      thisPerm=lapply(permutedTrt, '[', , b, drop=FALSE)
#      bcss=get.bcss(thisPerm)
      
      wmrpp.p=numeric(length(spar))
      for(s.i in seq_along(spar)){
#        wt=get.wt(bcss, spar[s.i])
        wdist=get.wdist(all.weights[s.i, b, ])
        tmp=.Call('mrppstats',wdist, permutedTrt, wtmethod, PACKAGE='MRPP')
        wmrpp.p[s.i]=mean(tmp[b]-tmp >=-eps)
      }
      plot(wmrpp.p~spar, main=b)
      stats[b]=min(wmrpp.p)
      if(b==1L){
        s0.i=max(which(spar==min(wmrpp.p)))
        s0=spar[s0.i]
        #wt0=get.wt(bcss, s0)
        wt0=all.weights[s0.i, 1L, ]
      }
    }
    
    #stats=.Call('mrppstats',y,permutedTrt, as.numeric(wtmethod[1L]), PACKAGE='MRPP')
    ans=list(statistic=c("Sparse Weighted MRPP statistic"=stats[1L]), all.statistics=stats, weights=wt0,
             p.value=mean(stats[1]-stats>=-eps), parameter=c("number of permutations"=B, 'group weight method'=wtmethod[1L], 'Smoothing'=s0, '#selected variables'=sum(wt0>0)),
             data.name=dname, #  .Random.seed=attr(permutedTrt,'.Random.seed'),  ## this is commented out since the random number seed in gmp:::urand.bigz will be used. 
             method=sprintf('%d-sample Sparse Weighted MRPP test',ntrt)
             )
    class(ans)='htest'
    ans
}

smrpp.test.formula <-
function(y, data, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    if(missing(data)) mf = model.frame(y) else mf <- model.frame(y, data=data)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    res=mf[[response]]
    ans=smrpp.test(res,trt=g, ...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " ~ "))
    if(!missing(data)) ans$data.name = paste(ans$data.name, "in data.frame", substitute(data))
    ans
}

smrpp.test <-
function(y,...) UseMethod("smrpp.test")
