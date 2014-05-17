if(FALSE){ ## this is the same a scaled t
dp7=function(x, center=0, a=1, m=1, log=FALSE)
{
    logans=-m*log(1+x*x/a/a)-log(a)-lbeta(m-.5,.5)
    if(isTRUE(log)) logans else exp(logans)
}
}
p7=function(x,y,method=c('both','PL','alpha'),lambdas=10^seq(1,13,length=30),alphas=10^(1:8), #lambda.n=10, 
    eps=.Machine$double.eps^.25, niter=1000L, verbose=TRUE)
{
    x=as.matrix(x); p=ncol(x); 
    ym=as.double(y-mean(y))
    xm=sweep(x,2,colMeans(x),'-')
    xnorm=sqrt(colSums(xm*xm))
    x.std=sweep(xm,2,xnorm,'/')
    beta0=drop(crossprod(x.std,ym))
    beta0=beta0*drop(crossprod(beta0)/crossprod(x.std%*%beta0))

    method=match.arg(method);
    if(method=='PL'){
        ans=.C('p7pmle',as.double(x.std), ym, length(ym), p, b=beta0,
				        alpha=1.0, as.double(eps),as.integer(verbose))
        ret=c(mean(y),ans$b/xnorm)
        attr(ret,'m')=.5*p/sum(1/(ans$alpha*ans$alpha/ans$b/ans$b+1))
        attr(ret,'alpha')=ans$alpha
        attr(ret,'sig2e')=mean((ym-x%*%ret[-1])^2)
        return(ret)
    }else if (method=='alpha'){
        ans=.C('p7alpha',as.double(x.std), ym, length(ym), p, b=rep(beta0,length(alphas)),
				        alpha=as.double((alphas)), length(alphas), 
                        as.double(eps),  as.integer(verbose))
        ret=rbind(mean(y), matrix(ans$b, p, length(alphas))/xnorm)
        attr(ret,'alpha')=ans$alpha
        return(ret)
    }

#    if(missing(lambdas)){
#        lambda.rg=range(beta0*beta0)/2
#        lambdas=10^seq(log10(lambda.rg[1]*.1), log10(lambda.rg[2]+.1), length=lambda.n)
#    }
    grids=expand.grid(sort(lambdas), rev(sort(alphas)))
    lambdas=grids[,1]
    alphas=grids[,2]
    lambda.n=length(lambdas)
#    print(cbind(log10(lambdas), log10(alphas)))


    ans=.C("p7", as.double(x.std), ym, length(ym), p, b=rep(beta0, length=p*lambda.n), #rep(beta0,length=p*lambda.n), 
                    as.double(lambdas), as.double(alphas), length(lambdas), as.double(eps),as.integer(niter),
                    as.integer(verbose))
    ret=rbind(mean(y),matrix(ans$b,p,length(lambdas))/xnorm)
    colnames(ret)=sprintf("lambda=%1.2E,alpha=%1.2E",lambdas,alphas)
    attr(ret,'alpha')=alphas
    attr(ret,'lambda')=lambdas
    resids=ym-x.std%*%ret[-1]
    attr(ret,'sig2e')=colMeans(resids*resids)
    attr(ret,'m')=lambdas/attr(ret,'sig2e')
    attr(ret,'df')=2*attr(ret,'m')-1
    attr(ret,'s')=alphas/sqrt(attr(ret,'df'))
    ret
}



if(FALSE){
    dyn.load("coordDesc.so")
    n=2000;
    p=9000;
#    p=10
    x=matrix(rnorm(n*p),n,p); 
    x=sweep(x,2,colMeans(x));
    x=sweep(x,2,sqrt(colSums(x*x)),'/');
    y=rnorm(n); 
    y=y-mean(y)
    beta0=crossprod(x,y)

    lambda=.5
    tmp7=p7(x,y,1,1,eps=1e-6)
    cor(y,x%*%tmp7[-1])
    cor(beta0,tmp7[-1])

    tmp7r=p7r(x,y,1,1)
    cor(y,x%*%tmp7r)
    cor(beta0,tmp7r)

    summary(tmp7[-1]-tmp7r)

    tmp7all=p7(x,y,lambda.n=3)
    tmp7all=p7(x,y,method='both', eps=1e-3)
    tmp7all=p7(x,y,method='alpha',eps=1e-3)
    tmp7all=p7(x,y,method='alpha',alphas=1,eps=1e-3)

    tmp7.a1=p7(x,y,method='both',alphas=1.0, eps=1e-3)
    
    tmp7.uni=apply(x,2,p7,y=y,method='both',lambdas=.1,alphas=.10,verbose=FALSE) ## discontinuous?
    plot(beta0,tmp7.uni[2,],pch='.',cex=2); abline(0,1,0,0,col=4)


    p7r=function(x,y,M=1, ALPHA=1, eps=.Machine$double.eps^.25, verbose=T) ## original slow R version for a single lambda
    {
        p=ncol(x)
        ym=y-mean(y)
        xm=sweep(x,2,colMeans(x),'-')
        xnorm=sqrt(colSums(xm*xm))
        x.std=sweep(xm,2,xnorm,'/')

    #    b=coef(lm.ridge(ym~x.std,lambda=.1))[-1];
        b=rep(0.0,p)
        iter=list(b)
        i=1
        repeat{
            last.b=b
            r=ym-x.std%*%b
            maxdiff=-Inf
            for(j in 1:p){
                r=r+x.std[,j]*b[j]
    #            lsb=drop(crossprod(x.std[,j],r))
#                cat(r,fill=T)
                lsb=sum(x.std[,j]*r);
                A=-lsb; B=2*M+ALPHA*ALPHA; C=-lsb*ALPHA*ALPHA
                Mterm=2*A*A*A-9*A*B+27*C; K=A*A-3*B; Nterm=(Mterm*Mterm-4*K*K*K)
                if(Nterm>0){
                    sqrtNterm=sqrt(Nterm);
                    b[j]=-1/3*(A+(abs(Mterm+sqrtNterm)/2)^(1/3)*ifelse(Mterm+sqrtNterm>0,1,-1)+
                             (abs(Mterm-sqrtNterm)/2)^(1/3)*ifelse(Mterm-sqrtNterm>0,1,-1))
                }else { ## all three roots are real
#                    theta=acos(-Mterm/2/sqrt(K*K*K))
#                    rts=2*sqrt(K)/3*cos((theta+(0:2)*2*pi)/3)+lsb/3
                    rts=Re(polyroot(c(C,B,A,1)))
                    cat(rts,fill=T)
                    obj=.5*(lsb-rts)^2+M*log(1+rts*rts/ALPHA/ALPHA)
                    b[j]=rts[which.min(obj)]
                }
                r=r-x.std[,j]*b[j]
                thisdiff=abs(b[j]-last.b[j])
                if(thisdiff>maxdiff) maxdiff=thisdiff
            }
            cat(sprintf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm))

            i=i+1
            iter[[i]]=b
    #        if(max(abs(last.b-b))<1e-6) break
            if(maxdiff<eps) break
            if(verbose)cat(cor(b,last.b),'\t',maxdiff,#max(abs(last.b-b)), 
    #                '\t',max(abs(last.b-b)/pmax(abs(last.b),abs(b)),na.rm=T),
                    fill=T)
        }
        final=b/xnorm
        attr(final,'iter')=iter
        final
    }



################   QTLMAS2010 dataset
load("/Users/longor/MRPP.influence/QTLMAS2010/raw.RData")

y=ydat[,2]; N=nrow(ydat)
X=t(apply(xdat[1:N,-1],1,function(xx)as.integer(colSums(matrix(xx,2))-3)))  ## design mat

noninformative.snp=which(sd(X)==0)
X=X[,-noninformative.snp]
dim(X)
#[1] 2326 9768

#y.bar=mean(ydat[,2])
#X.bar=colMeans(X)
#X.norm=sqrt(colSums(sweep(X,2,X.bar,'-')^2))
#X.std=sweep(sweep(X,2,X.bar,'-'),2,X.norm,'/')[,-noninformative.snp]
#dim(X.std)

p7.fit=p7(X, y)

}