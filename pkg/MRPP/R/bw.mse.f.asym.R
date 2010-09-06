bw.mse.F.asym <-
function(x,iter=TRUE,eps=1e-3,iter.max=100,verbose=FALSE) ## cubic
{require(ks)
        K1=0.2820947918
        bw0=bw.nrd(x)
        

        n.iter=1
        repeat{
            f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
            F.1=mean(pnorm(x[1]-x,0,bw0))
            df=drvkde(x,1,bw0,se=FALSE)
            df.1=approx(df$x.grid[[1]],df$est,xout=x[1])$y
            bw=f.1*K1/df.1/(.5-F.1)
            if(!iter || abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
            if(verbose)cat(bw,fill=T)
            bw0=bw
            n.iter=n.iter+1
        }
}

