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

if(FALSE){
mrpp.test.default <-
function(y, ...) {
    ans=mrpp.test.dist(dist(y),...)
    repl.text=paste("Response data ", deparse(substitute(y)), sep='')
    ans$data.name=gsub("\"dist\" object dist(y)", repl.text, ans$data.name, fixed=TRUE)
    ans
}

mrpp.test.formula <-
function(y, data, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    if(missing(data)) mf = model.frame(y) else mf <- model.frame(y, data=data)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    res=mf[[response]]
    ans=mrpp.test(dist(res),trt=g, ...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " ~ "))
    if(!missing(data)) ans$data.name = paste(ans$data.name, "in data.frame", substitute(data))
    ans
}
}

smrpp.test <-
function(y,...) UseMethod("smrpp.test")

