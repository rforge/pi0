bw.mse.pdf.asym=function(x,iter=FALSE,eps=1e-3,iter.max=100,verbose=FALSE)
{require(ks)

        bw0=bw.nrd(x)
        Rkern=density(x,bw0,from=x[1],to=x[1],n=1,give.Rkern=TRUE)

        n.iter=1
        repeat{
            f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
            ddf=drvkde(x,2,bw0,se=FALSE)
            ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
            bw=(f.1/ddf.1/ddf.1*Rkern/length(x))^.2
            if(!iter || abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
            if(verbose)cat(bw,fill=TRUE)
            bw0=bw
            n.iter=n.iter+1
        }
}
