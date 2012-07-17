bw.mse.pdf.asym=function(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, verbose=FALSE)
{#require(ks)
    if(is.function(start.bw)) bw0=start.bw(x)
    if(is.numeric(start.bw)) bw0=start.bw
    if(is.character(start.bw)) bw0=call(start.bw, x)

    Rkern=density(x,bw0,from=x[1L],to=x[1L],n=1L,give.Rkern=TRUE)

    n.iter=1L
    xdiff=x[1L]-x
    repeat{
#### these lines are using the ks:::drvkde and or ks::kdde functions and  are replaced by explicit calculations
#        f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
#        ddf=drvkde(x,2,bw0,se=FALSE)
#        ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
            tmp=xdiff/bw0
            f.1=mean(dnorm(tmp))/bw0
            ddf.1=mean(dnorm(tmp)*(tmp*tmp-1))/bw0/bw0/bw0
        bw=(f.1/ddf.1/ddf.1*Rkern/length(x))^.2
        if(abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
        if(verbose && isTRUE(n.iter%%verbose==0))cat(bw,fill=TRUE)
        bw0=bw
        n.iter=n.iter+1
    }
}

if(FALSE){
    f.1=kdde(x, bw0, deriv.order=0, eval.points=x[1])$estimate
    ddf.1=kdde(x, bw0, deriv.order=2, eval.points=x[1])$estimate

    locpoly(x,  drv = 0L, degree=0, kernel = "normal", gridsize=2L, range.x=0:1+x[1], bandwidth=bw0)$y

    f.1=mean(dnorm((x[1]-x)/bw0))/bw0


        f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
        ddf=drvkde(x,2,bw0,se=FALSE)
        ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
        f.1;ddf.1
}
