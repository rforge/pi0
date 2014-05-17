#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>

#define negLarge (-DOUBLE_XMAX)    
#define posLarge (DOUBLE_XMAX)    

#define oneThird (0.3333333333333333333)
#define Factor (1.9802913004322129507587985355868732882822056101876)
#define log1pff (1.5936242600400400923230418758751602417890024248188)
#define ffp1df (2.48526751266005208052) 
// (row, col)th element of matrix mat that has nrow rows, indexing from 0. 
#define mat0(mat, row, col, nrow) (*((mat)+(row)+(col)*(nrow)))
#define cubeRoot(x) ((sign(x)) * (R_pow(abs(x), oneThird)))

double R_INLINE lassotFit1(double betahat, double lambda, double alpha)
{
	double root,  A, B, C, M, sqrtN, K;
	root=fabs(betahat)-lambda; 
	root=sign(betahat)* (root>0.0 ? root : 0.0);
	if(fabs(root) <= Factor*alpha) return root;
	
	A=-betahat; B=alpha*alpha + lambda*alpha * ffp1df; C=-alpha*alpha*betahat;
	K=A*A-3.0*B; M=2.0*A*A*A-9.0*A*B+27.0*C; sqrtN=sqrt(M*M-4.0*K*K*K);
	// note: this is assuming N>0 
	return -oneThird * (A+cubeRoot(.5*(M+sqrtN))+cubeRoot(.5*(M-sqrtN)));  
}

void lassot(double *xstd, double *ym, int *n, int *p, double *beta, 
				double *lambda, double *alpha, int *nlambda, double *eps, int *niter, int *verbose)
// double exponential log prior for |x|<=Factor*alpha; BayesA type log prior for |x|>Factor*alpha
// Assuming continuity condition (convexness) met and alpha and lambda all being positive
// upon return, lambda will contain sum of coordinate-wise degrees of freedom
{	
	double  maxdiff;
	double  *resid, *dfs;
	int j, ilambda, iter, tmpint, tmpint1;
	double A, B, C, K, Mterm,  sqrtNterm, tmpp, tmpq,
			tmpDouble, tmpDouble1, tmpDouble2,
			alphaFactor, lastLambda, lastAlpha;
	char cN='n';
	double negOne=-1.0;
	double done=1.0;
	int ione=1;

	// initialization
	lastLambda = -1.0; lastAlpha = -1.0;
	resid   = (double*) Calloc(*n, double);
	dfs     = (double*) Calloc(*nlambda, double);
	
	for(ilambda=0; ilambda<*nlambda; ++ilambda){
		// warn start beta
		if(*(lambda+ilambda) == lastLambda ){
			memcpy(beta+(*p)*ilambda, beta+(*p)*(ilambda-1), sizeof(double)*(*p));
			if( *(alpha+ilambda) == lastAlpha ) continue;
		}else{
			lastLambda = *(lambda+ilambda); 
			if(ilambda){
				for(tmpDouble=posLarge, tmpint=0; tmpint<ilambda; ++tmpint) 
					if(fabs(*(alpha+tmpint) - *(alpha+ilambda)) < tmpDouble) tmpDouble = fabs(*(alpha+tmpint) - *(alpha+ilambda)) ;
				for(tmpDouble1=posLarge, tmpint1=-1, tmpint=0; tmpint<ilambda; ++tmpint)
					if(fabs(*(alpha+tmpint) - *(alpha+ilambda)) == tmpDouble &&
					   fabs(*(lambda+tmpint) - *(lambda+ilambda)) < tmpDouble1)
						tmpint1=tmpint, tmpDouble1 = fabs(*(lambda+tmpint) - *(lambda+ilambda)); 
				if(tmpint1 < 0) tmpint1 = ilambda-1L;
				memcpy(beta+(*p)*ilambda, beta+(*p)*tmpint1, sizeof(double)*(*p));
			}
		}
		
		lastAlpha = *(alpha+ilambda);
		
		// init resid
		memcpy(resid, ym, sizeof(double)*(*n));
		F77_CALL(dgemv) ( &cN, n, p, &negOne, xstd, n, beta+(*p)*ilambda, &ione, &done, resid, &ione ); // resid=ym-xstd%*%beta[,ilambda]

		// iterate
		maxdiff=0.0;
		for(iter=0; iter<*niter; ++iter){
			for(j=0; j<*p; ++j){
				tmpDouble = lassotFit1(F77_CALL(ddot)(n, xstd+(*n)*j, &ione, resid, &ione) + mat0(beta, j, ilambda, *p), lastLambda, lastAlpha); 
				tmpDouble1 = mat0(beta, j, ilambda, *p) - tmpDouble ;
				if( tmpDouble1 != 0.0 ){
					F77_CALL(daxpy)(n, &tmpDouble1, xstd+(*n)*j, &ione, resid, &ione);
					if(fabs(tmpDouble1) > maxdiff) maxdiff = fabs(tmpDouble1);
					mat0(beta, j, ilambda, *p) = tmpDouble; 
				}
			}// of j (variable)
			if(*verbose)	Rprintf("lambda=%.3f\talpha=%.3f\titer=%d\t|max.diff|=%.15f\n",lastLambda, lastAlpha, iter, maxdiff);
			
			if(maxdiff <= *eps) break;
			maxdiff=0.0;
			R_CheckUserInterrupt();
		} // of iter
		
		// computing coordinate-wise degrees of freedom
		dfs[ilambda]=0.0;
		alphaFactor = Factor * lastAlpha; 
		for(j=0; j<*p; ++j){
			tmpDouble=mat0(beta, j, ilambda, *p); // beta[j,i]
			if(tmpDouble == 0.0) continue;
			if(tmpDouble <= alphaFactor) {dfs[ilambda]+=1.0;  continue;}
			tmpDouble2 = F77_CALL(ddot)(n, xstd+(*n)*j, &ione, resid, &ione) + tmpDouble; //  bhat
			A = -tmpDouble2; B=lastAlpha*lastAlpha+lastLambda*lastAlpha*ffp1df; C=-lastAlpha*lastAlpha*tmpDouble2;
			tmpp=2.0*tmpDouble2*tmpDouble2 - 3.0*lastAlpha*lastLambda*ffp1df + 6.0*lastAlpha*lastAlpha;
			tmpq=8.0*tmpDouble2*tmpDouble2 -R_pow(ffp1df*lastLambda,2.0)-20.0*lastAlpha*lastLambda*ffp1df+8.0*lastAlpha*lastAlpha;
			Mterm=2.0*R_pow(A,3.0)-9.0*A*B+27.0*C ;
			K=A*A-3.0*B; sqrtNterm=sqrt(Mterm*Mterm -4.0*R_pow(K,3.0));
			tmpDouble1=9.0*lastAlpha*lastAlpha*tmpDouble2*tmpq/sqrtNterm;
			dfs[ilambda]+=((tmpp-tmpDouble1)/cubeRoot(2.0*R_pow(Mterm+sqrtNterm,2.0)) + 
						   (tmpp+tmpDouble1)/cubeRoot(2.0*R_pow(Mterm-sqrtNterm,2.0)) + 
						   1.0) / 3.0;
		}
	} // ilambda
	
	memcpy(lambda, dfs, sizeof(double)*(*nlambda));
	Free(dfs);
	Free(resid);
}

/* ///  OLD CODE
	thislambdaN=0;
	for(i=0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=mat0(x, i, j, *n)*mat0(beta, j, 0, *nlambda);
//		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=mat0(x, i, (*p)-1, *n)*mat0(beta, (*p)-1, 0, *nlambda);


//	lastaa=0.0;
	lastlambda=-1;
	for(ilambda=(*nlambda)-1; ilambda>=0; --ilambda) {
		thislambdaN=ilambda*(*p);
		aa=alpha[ilambda]*alpha[ilambda];
//		if(aa==lastaa)
//			memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nlambda)), sizeof(double)*(*p)); // hotstart
//		lastaa=aa;
		if(lambda[ilambda]==lastlambda)
			memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nlambda)), sizeof(double)*(*p)); // hotstart
		lastlambda=lambda[ilambda];

		alphaFactor=alpha[ilambda]*Factor;
		lambdaprime=lambda[ilambda]*log1pff/alphaFactor;
		iter=0;
		do{
			maxdiff=maxreldiff=0.0;  //negLarge; //nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
				lastbeta=beta[thislambdaN+j];
				tmproot=sign(lsb)*(fabs(lsb)>lambdaprime ? fabs(lsb)-lambdaprime : 0.0);
				if(abs(tmproot)>=alphaFactor){
					A=-lsb; B=2.0*lambda[ilambda]+aa; C=-lsb*aa;
					Mterm=2.0*A*A*A-9.0*A*B+27.0*C; K=A*A-3.0*B; Nterm=(Mterm*Mterm-4.0*K*K*K);
					if(Nterm>0.0){
						sqrtNterm=sqrt(Nterm);
						tmproot=-oneThird*(A+R_pow(fabs(Mterm+sqrtNterm)*.5, oneThird)*sign(Mterm+sqrtNterm)+
										 R_pow(fabs(Mterm-sqrtNterm)*.5, oneThird)*sign(Mterm-sqrtNterm));
	//					tmproot=fabs(tmproot)>=alphaFactor ? tmproot : DOUBLE_XMAX;
					}else{  // three real roots
						sqrtK=sqrt(K);
						theta=acos(-.5*Mterm/sqrtK/K);
						root1=2.0*oneThird*sqrtK*cos(theta*oneThird)+lsb*oneThird;
	//					root1=fabs(root1)>=alphaFactor ? root1:DOUBLE_XMAX;
						
						root2=2.0*oneThird*sqrtK*cos((M_2PI+theta)*oneThird)+lsb*oneThird;
	//					root2=fabs(root2)>=alphaFactor	? root2:DOUBLE_XMAX;

						root3=2.0*oneThird*sqrtK*cos((2.0*M_2PI+theta)*oneThird)+lsb*oneThird;
	//					root3=fabs(root3)>=alphaFactor	? root3:DOUBLE_XMAX;

						root1=.5*(lsb-root1)*(lsb-root1)+lambda[ilambda]*log1p(root1*root1/aa) <
							  .5*(lsb-root2)*(lsb-root2)+lambda[ilambda]*log1p(root2*root2/aa) ?
							  root1 : root2 ;
						tmproot=
							 .5*(lsb-root1)*(lsb-root1)+lambda[ilambda]*log1p(root1*root1/aa) <
							 .5*(lsb-root3)*(lsb-root3)+lambda[ilambda]*log1p(root3*root3/aa) ?
							 root1 : root3 ;
					}
					if(fabs(tmproot)<alphaFactor) tmproot=DOUBLE_XMAX;
				}
				beta[thislambdaN+j]=tmproot;
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff)
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 
				if(beta[thislambdaN+j]==0.0){
					if(lastbeta!=0.0 && maxreldiff<1.0) maxreldiff=1.0;
				}else
					if(fabs(beta[thislambdaN+j]-lastbeta)/fabs(beta[thislambdaN+j])>maxreldiff) 
					maxreldiff=fabs(beta[thislambdaN+j]-lastbeta)/fabs(beta[thislambdaN+j]) ; 
			}
			if(*verbose)printf("lambda=%4.2E\talpha=%4.2E\tmaxdiff=%1.15f\tmaxreldiff=%1.15f\n",lambda[ilambda],alpha[ilambda],maxdiff,maxreldiff);
			R_CheckUserInterrupt();
		}while((maxdiff>*eps || maxreldiff>0.1) && ++iter<*niter);
	}
	free(resid);
}

*/
