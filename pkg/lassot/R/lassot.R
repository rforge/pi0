lassotFactor=1.98029130043221295 #remaining digits are 07587985355868732882822056101877  # explicit sol: lassotFactor=sqrt(exp(2+LambertW::W(-2/exp(1)^2))-1)
ffp1=1+lassotFactor^2
ffp1df=ffp1/lassotFactor
log1pff=log(ffp1)
continuityFactor=(lassotFactor^2-1)/lassotFactor/(lassotFactor^2+1)
log2=log(2)
oneThird=1/3
cubeRoot=function(x)sign(x)*abs(x)^oneThird

dlassot1=function(x,a,m,df=2*m-1,s=a/sqrt(df), log=FALSE) {  ## only when lambda=1
    half.const=pt(-lassotFactor*sqrt(df),df)*a*beta(df/2,.5)+a*lassotFactor/m/log1pff*(1-(1+lassotFactor*lassotFactor)^-m)
    logans=ifelse(abs(x)<=lassotFactor*a,
        -abs(x)*m*log1pff/lassotFactor/a, 
        -m*log(1+x*x/a/a))-log(half.const)-log(2)
    if(isTRUE(log))logans else exp(logans)
}

dlassot=function(x, rate=1/continuityFactor/ffp1df, alpha=2/ffp1df, log=FALSE)
{
    a=alpha
	stopifnot(all(rate>0, alpha>1/ffp1df))
	
	# print((df=a*ffp1df-1))
	# print((s=a/sqrt(df)))
    #a=(1+df)/ffp1df
	#m=(1+df)/2
	lambda=rate
	f=lassotFactor
	const=1/((2 - 2 *exp(-a *f* lambda))/lambda + 1/(-f + a *(1 + f^2)* lambda)*2* a* f^(2 - (a* (1 + f^2) *lambda)/f)     *hypergeo::hypergeo(  (a *(1 + f^2) *lambda)/(2* f),    1/2* (-1 + (a* lambda)/f + a* f* lambda), (f + a* (1 + f^2) *lambda)/(  2 *f), -(1/f^2)   ) )
	
	logans=log(Re(const))+ ifelse(abs(x)<=lassotFactor*a,
        -abs(x)*rate,  
        -ffp1df/2*rate*a*log(1+x*x/a/a))
	if(isTRUE(log)) logans else exp(logans)
}

lassot.fit=function(x,y,lambdas=10^seq(1,5,length=10),alphas=1.5^seq(1,17,length=10), 
    method=c('Coordinate','Newton','R'), eps=.Machine$double.eps^.25*ncol(x), niter=1000L, verbose=TRUE)
{
	continuity=TRUE  ## FALSE case not fully implemented
	
    x=as.matrix(x); p=ncol(x);     n=length(y)
    ym=as.double(y-mean(y))/sd(y)
    xm=sweep(x,2L,colMeans(x),'-')
    xnorm=sqrt(colSums(xm*xm))
    x.std=sweep(xm,2L,xnorm,'/')
    beta0=drop(crossprod(x.std,ym))
    beta0=beta0*drop(crossprod(beta0)/crossprod(x.std%*%beta0))

    grids=expand.grid(rev(sort(alphas)),rev(sort(lambdas)) )
    lambdas=grids[,2]
    alphas=grids[,1]
	if(isTRUE(continuity)) alphas=pmax(continuityFactor*lambdas, alphas)
    lambda.n=length(lambdas)

    method=match.arg(method)
    if(method=='Coordinate'){
        ans=.C("lassot", as.double(x.std), ym, length(ym), p, 
                    b=rep(0.0, length=p*lambda.n), #rep(beta0,length=p*lambda.n), 
                    lambda=as.double(lambdas), as.double(alphas), length(lambdas), as.double(eps),
                    as.integer(niter), as.integer(verbose))
		attr(ans$b, 'df')=ans$lambda
    }else if(method=='R') {
		ans=lassot.fitR(x.std, ym, lambdas, alphas, continuity, eps, niter,verbose)
		ans=list(b=ans)
	}else if (method=='Newton') { ## too slow
        xpx=crossprod(x.std)
        obj=function(betas, lambda, alpha)
            .5*sum(crossprod(ym-x.std%*%betas))+
            lambda*sum(ifelse(abs(betas)<=alpha, log(1+betas*betas/alpha/alpha), abs(betas)/a-1+log(2)))
        der=function(betas, lambda, alpha)
            -beta0+xpx%*%betas+lambda*ifelse(abs(betas)<=alpha,2*betas/(betas*betas+alpha*alpha),
                            ifelse(betas>0,1,-1)/alpha)
        ret=matrix(mean(y), p+1, lambda.n)
        for(i in 1:lambda.n)
            ret[-1,i]=nlminb(rep(0,p),obj,der,lambda=lambdas[i],alpha=alphas[i],control=list(trace=1))$par
    }
#    colnames(ret)=sprintf("lambda=%1.2E,alpha=%1.2E",lambdas,alphas)
    # plogLik=function(betas, a, lambda){sig2=crossprod(ym-x.std%*%betas)/n
        # m=lambda/sig2
        # sum(dnorm(ym,x.std%*%betas,sqrt(sig2),log=TRUE))+sum(dlassot(betas,lambda,a,log=TRUE))
    # }
    betas=matrix(ans$b,p,lambda.n)
    # pll=numeric(lambda.n)
    # for(i in 1:lambda.n)pll[i]=plogLik(betas[,i],alphas[i],lambdas[i])
    # ibest=which.max(pll) # this never return NaN/NA :)
	
	sse=colSums((ym-x.std%*%betas)^2)
	gcv=sse/n/(1-attr(ans$b, 'df')/n)^2
	ibest=which.min(gcv)
    ret=c(mean(y),betas[,ibest]/xnorm*sd(y))
    attr(ret,'alpha')=alphas[ibest]
    attr(ret,'lambda')=lambdas[ibest]
    attr(ret,'sig2e')=crossprod(y-x%*%ret[-1]-mean(y))/(n-min(n-1,attr(ans$b,'df')[ibest]))
    attr(ret, 'df')=attr(ans$b,'df')[ibest]
    # attr(ret,'dfs')=attr(ans$b, 'dfs')
    # attr(ret,'t.df')=alphas[ibest]*ffp1df-1
    # attr(ret,'s')=alphas[ibest]/sqrt(attr(ret,'t.df'))
    attr(ret,'gcv')=gcv
    attr(ret,'betas')=betas/xnorm*sd(y)
    attr(ret,'dfs')=attr(ans$b, 'df')
	attr(ret,'lambdas')=lambdas
    attr(ret,'alphas')=alphas
    ret
}

lassot.fit1R=function(beta.hat, lambda, alpha, df=FALSE)
{   
    pls=function(beta).5*(beta-beta.hat)^2+lambda*ifelse(abs(beta)<=lassotFactor*alpha, abs(beta), ffp1df/2*alpha*log(1+beta^2/alpha^2))
	root0=sign(beta.hat)*pmax(abs(beta.hat)-lambda, 0)
	if(abs(beta.hat)<=lambda+lassotFactor*alpha) {
		if(isTRUE(df))	attr(root0, 'df')=as.numeric(root0!=0)
		return(root0)
	}
	A=-beta.hat
	B=alpha^2+lambda*alpha*ffp1df
	C=-alpha^2*beta.hat
	if(FALSE){ ## this part is correct but using polyroot
		sols=polyroot(c(C,B,A,1))
		sol=Re(sols[abs(Im(sols))<1e-6])
		obj=pls(sol)
		ans=sol[which.min(obj)]
		return(ans)
	}
	
	K=A*A-3*B
	M=2*A^3-9*A*B+27*C  ## = 9*alpha*(lambda*ffp1df-2*alpha)*beta.hat-2*beta.hat^3
	N=M*M-4*K^3; sqrtN=sqrt(N)
	if(N>0) {
		ans=-oneThird*(A+cubRoot((M+sqrtN)/2)+cubRoot((M-sqrtN)/2))
		if(isTRUE(df)){
			tmpp=2*beta.hat^2 - 3*alpha*lambda*ffp1df + 6*alpha^2
			tmpq=8*beta.hat^2 - (ffp1df*lambda)^2 - 20*alpha*lambda*ffp1df + 8*alpha^2
			ninea2bqdsqrtN=9*alpha^2*beta.hat*tmpq/sqrtN
			attr(ans,'df')=((tmpp-ninea2bqdsqrtN)/cubeRoot(2*(M+sqrtN)^2) + (tmpp+ninea2bqdsqrtN)/cubeRoot(2*(M-sqrtN)^2) + 1)/3
		}
		return (  ans  )
	}
	
	browser() ## this line should not be reached when continuity requirement is met, because the obj is convex
	theta=acos(-M/2/K/sqrt(K))
	ell=0:2
    insideTrig=(theta+2*ell*pi)/3
	sol=2/3*sqrt(K)*cos(insideTrig)+oneThird*beta.hat
	if(FALSE){
		dKdBetahat=2*beta.hat
		dThetadBetahat=-3*alpha*ffp1df*beta.hat/4/K^2.5/sqrt(1-(alpha*ffp1df*lambda)^2/16/K^3)
		dSoldBetahat=(3*cos(insideTrig)*dKdBetahat-2*K*dThetadBetahat*sin(insideTrig)+3*sqrt(K))/9/sqrt(K)
	}
	obj=pls(sol)
	ans=sol[which.min(obj)]
	attr(ans, 'obj')=obj
#	attr(ans, 'deriv')=dSoldBetahat
	ans
}

lassot.fitR=function(x.std, ym, lambdas, alphas, continuity=TRUE, eps=.Machine$double.eps^.25*ncol(x.std), niter=1000L,verbose=TRUE)
{
	if(isTRUE(continuity)) alphas=pmax(alphas, continuityFactor*lambdas)
	nlambdas=length(lambdas)
	n=length(ym)
	p=ncol(x.std)
	betas=dfs=matrix(0, p, nlambdas)
	# selected=c()
	# nselected=length(selected)
	# xpy=crossprod(x.std, ym)
	# xpx=matrix(numeric(0L), p, 0L)
	lastLambda=-Inf
	lastAlpha=-Inf
	for(i in seq(nlambdas)) {
		### warm start
		if(lambdas[i]==lastLambda) {
			betas[,i]=betas[,i-1L]
			if( alphas[i]==lastAlpha) next
		}else {
			lastLambda=lambdas[i]
			minAlphaDiff=suppressWarnings(min(abs(alphas[varComp::safeseq(i-1L)]-alphas[i])))
			tmp.idx=which(abs(alphas[varComp::safeseq(i-1L)]-alphas[i])==minAlphaDiff)
			tmp.i=which.min(abs(lambdas[tmp.idx]-lambdas[i]))
			if(length(tmp.i)>0L) {
				betas[,i]=betas[,tmp.idx[tmp.i]]
			}else if(i>1L) betas[,i]=betas[,i-1L] 
		}
		
		lastAlpha=alphas[i]
		
		resids=ym-x.std%*%betas[,i]
		betas.old=betas[,i]
		for(iter in seq(niter)){
			for(j in seq(p)){
				if(TRUE){
					ans=lassot.fit1R(sum(x.std[,j]*resids)+betas[j,i], lambdas[i], alphas[i])
					if(ans!=betas[j,i]) {
						resids=resids-x.std[,j]*(ans-betas[j,i])
						betas[j,i]=ans
					}
				}
				
				if(FALSE){ ## this also works; slightly slower
					if(betas[j,i]!=0) resids=resids + x.std[,j]*betas[j,i]
					betas[j,i]=lassot.fit1R(sum(x.std[,j]*resids), lambdas[i], alphas[i])
					if(betas[j,i]!=0) resids=resids-x.std[,j]*betas[j,i]
				}
				
				if(FALSE){ ## the so-called covariance method, not tested yet; not expected to be working well for large p
					beta.hat=xpy[j]-crossprod(xpx[j,], betas[selected])+betas[j,i]
					ans=lassot.fit1R(beta.hat, lambdas[i], alphas[i])
					if(ans!=0 && betas[j,i]==0){
						selected=sort(c(selected, j))
						xpx=crossprod(x.std, x.std[,selected,drop=FALSE])
					}
					betas[j,i]=ans
				}
			}
			mdiff=max(abs(betas.old-betas[,i]))
			if(isTRUE(verbose)) cat('lambda=',lambdas[i],'\talpha=',alphas[i],'\titer=', iter, '\tmax.diff=', mdiff, '\n')
			if(mdiff<=eps) break
			betas.old=betas[,i]
		}
		## computing coordinate-wise degrees of freedom
		lasso.idx=(abs(betas[,i])<=lassotFactor*alphas[i])
		dfs[lasso.idx,i]=as.numeric(betas[lasso.idx,i]!=0)
		for(j in which(!lasso.idx)){
			bhat=sum(x.std[,j]*resids)+betas[j,i]
			A=-bhat; B=alphas[i]^2+lambdas[i]*alphas[i]*ffp1df; C=-alphas[i]^2*bhat
			tmpp=2*bhat^2 - 3*alphas[i]*lambdas[i]*ffp1df + 6*alphas[i]^2
			tmpq=8*bhat^2 - (ffp1df*lambdas[i])^2 - 20*alphas[i]*lambdas[i]*ffp1df + 8*alphas[i]^2
			M=2*A^3 - 9*A*B + 27*C; K = A^2 - 3*B; N = M^2 - 4 *K^3; sqrtN=sqrt(N)
			ninea2bqdsqrtN=9*alphas[i]^2*bhat*tmpq/sqrtN
			dfs[j,i]=((tmpp-ninea2bqdsqrtN)/cubeRoot(2*(M+sqrtN)^2) + (tmpp+ninea2bqdsqrtN)/cubeRoot(2*(M-sqrtN)^2) + 1)/3
		}
		
	}
	ans=betas
	# attr(ans, 'dfs')=dfs
	attr(ans, 'df')=colSums(dfs)
	ans
}
					
if(FALSE) {
    dyn.load("coordDesc.so")
	
	library(glmnet); library(lars)
    n=200;
    p=900;
#    p=10
	set.seed(3942)
    x=matrix(rnorm(n*p),n,p); 
    x=sweep(x,2,colMeans(x));
    x=sweep(x,2,sqrt(colSums(x*x)),'/');
    y=rnorm(n,0,2)+x%*%c(1:5, rep(0,p-5)); 
    y=y-mean(y)
    beta0=crossprod(x,y)
	lambda.max=max(abs(beta0))

    # lambda=.5
    # tmp7=p7(x,y,1,1,eps=1e-6)
    # cor(y,x%*%tmp7[-1])
    # cor(beta0,tmp7[-1])

    # tmp7r=p7r(x,y,1,1)
    # cor(y,x%*%tmp7r)
    # cor(beta0,tmp7r)

    # summary(tmp7[-1]-tmp7r)

    # tmp7all=p7(x,y,lambda.n=3)
    # tmp7all=p7(x,y,method='both', eps=1e-3)
    # tmp7all=p7(x,y,method='alpha',eps=1e-3)
    # tmp7all=p7(x,y,method='alpha',alphas=1,eps=1e-3)

    # tmp7.a1=p7(x,y,method='both',alphas=1.0, eps=1e-3)
    
    # tmp7.uni=apply(x,2,p7,y=y,method='both',lambdas=.1,alphas=.10,verbose=FALSE) ## discontinuous?
    # plot(beta0,tmp7.uni[2,],pch='.',cex=2); abline(0,1,0,0,col=4)

	lasso.fit=glmnet(x,y); tmp.df=sum(0!=lasso.fit$beta[,25]); lasso.bet=lasso.fit$beta[,25]
	lars.fit=lars(x,y, use.Gram=FALSE); lars.bet=lars.fit$beta[tmp.df+1,]
    lassot1=lassot.fit(x,y, lambdas=10^seq(log10(lambda.max), log10(lambda.max)-2, length=100), alphas=1000,method='R', eps=1e-7*p); 
		lassot1.bet=attr(lassot1, 'betas')[,which.min(abs(attr(lassot1, 'dfs')-tmp.df))]
	plot(beta0, lasso.bet)
		points(beta0, lars.bet, col=2, pch=2)
		points(beta0, lassot1.bet, col=4, pch=3)
    lassot1c=lassot.fit(x,y, lambdas=10^seq(log10(lambda.max), log10(lambda.max)-2, length=100), alphas=1000,eps=1e-7*p); 
		lassot1c.bet=attr(lassot1c, 'betas')[,which.min(abs(attr(lassot1c, 'dfs')-tmp.df))]
		points(beta0, lassot1c.bet, col=5, pch=5)
	
	
    lassotall=lassot.fit(x[,1],y,method='R')
    lassotall=lassot.fit(x,y, lambdas=10^seq(-2,1,length=10),alphas=1.5^seq(1,17,length=20),method='R')
    lassot.uni=apply(x,2,lassot,y=y,method='Coordinate',lambdas=1.0,alphas=1.0,verbose=F) ## discontinuous?
    plot(beta0,lassot.uni,pch='.',cex=2);abline(0,1,0,0,col=4)
    lassot.uni=apply(x,2,lassot,y=y,method='Coordinate',lambdas=.10,alphas=.10,verbose=F) ## discontinuous?
    plot(beta0,lassot.uni,pch='.',cex=2);abline(0,1,0,0,col=4)



################   QTLMAS2010 dataset
load("/Users/longor/MRPP.influence/QTLMAS2010/raw.RData")

Y=ydat[,2]; N=nrow(ydat)
X=t(apply(xdat[1:N,-1],1,function(xx)as.integer(colSums(matrix(xx,2))-3)))  ## design mat

noninformative.snp=which(sd(X)==0)
X=X[,-noninformative.snp]
dim(X)
#[1] 2326 9768
pred.X=t(apply(xdat[-(1:N),-1],1,function(xx)as.integer(colSums(matrix(xx,2))-3)))  
pred.X=pred.X[,-noninformative.snp]

#y.bar=mean(ydat[,2])
#X.bar=colMeans(X)
#X.norm=sqrt(colSums(sweep(X,2,X.bar,'-')^2))
#X.std=sweep(sweep(X,2,X.bar,'-'),2,X.norm,'/')[,-noninformative.snp]
#dim(X.std)
beta.uni=drop(crossprod(X,Y))/colSums(X*X)
beta.uni=beta.uni*drop(crossprod(X%*%beta.uni,Y)/crossprod(X%*%beta.uni))
sig2e=crossprod(Y-X%*%beta.uni)/N
df=4

lassot.fit=lassot(X, Y, lambdas=10^seq(3.5,4.5,length=10), alphas=1.5^seq(5,20,length=50),eps=1e-6,niter=5e3)

cor(truth[1:nrow(pred.X)+nrow(ydat),'NA.QT'],pred.X%*%lassot.fit[-1])
#          [,1]
#[1,] 0.6357595

cors=cor(truth[1:nrow(pred.X)+nrow(ydat),'NA.QT'],pred.X%*%attr(lassot.fit,'betas'))
tail(sort(cors))
#[1] 0.5454816 0.5474577 0.5477132 0.5498540 0.5689820 0.5741520

plot(log10(attr(lassot.fit,'lambdas')), cors)
plot(log(attr(lassot.fit,'alphas'),1.5), cors)

lassot.lam=lassot(X, Y, lambdas=10^seq(2.5,5.5,length=50), alphas=50*4,eps=1e-6,niter=5e3)
tail(sort(cor(truth[1:nrow(pred.X)+nrow(ydat),'NA.QT'],pred.X%*%attr(lassot.lam,'betas'))))

}
