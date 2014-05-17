// unnecessary for release
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<unistd.h>

// necessary
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#define negLarge (-DOUBLE_XMAX)    

void logabs(double *x, double *y, int *n, int *p, double *beta, 
				double *lambda, int *nlambda, double *eps, int *verbose)
// log(abs(beta)) penalty
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ilambda, thislambdaN, nzero;
	double lsb, delta;

	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=0;
	for(i=0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
//		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];


	for(ilambda=(*nlambda)-1; ilambda>=0; --ilambda) {
		thislambdaN=ilambda*(*p);
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nlambda)), sizeof(double)*(*p)); // hotstart
		do{
			maxdiff=negLarge; nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
				/** this is the only block that may change for different priors **/
				delta=lsb*lsb-2.0*lambda[ilambda];
				if(lastbeta=beta[thislambdaN+j])nzero++;
				beta[thislambdaN+j]=delta>0.0? .5*(lsb+sign(lsb)*sqrt(delta)):0.0;
				/** end: possible changes **/
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 
			}
			if(verbose)printf("lambda=%1.2E\tmaxdiff=%1.15f\n",lambda[ilambda],maxdiff);
			R_CheckUserInterrupt();
		}while(maxdiff*nzero>*eps);
	}
	free(resid);
}


void logabsp1(double *x, double *y, int *n, int *p, double *beta, 
				double *lambda, int *nlambda, double *eps, int *verbose)
// log(abs(beta)+1) penalty
{
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ilambda, thislambdaN, nzero;
	double lsb, cutoff;

	resid   =(double *)malloc(sizeof(double)*(*n));
	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=0;
	for(i=0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
//		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];


	for(ilambda=(*nlambda)-1; ilambda>=0; --ilambda) {
		thislambdaN=ilambda*(*p);
		cutoff=lambda[ilambda]<=1.0 ? lambda[ilambda] : 2.0*sqrt(lambda[ilambda])-1.0;
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nlambda)), sizeof(double)*(*p)); // hotstart
		do{
			maxdiff=negLarge; nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
				if(lastbeta=beta[thislambdaN+j])nzero++;
				beta[thislambdaN+j]=fabs(lsb)>cutoff ? 
					.5*(lsb+sign(lsb)*(sqrt((lsb+sign(lsb))*(lsb+sign(lsb))-4.0*lambda[ilambda])-1.0)) : 0.0;
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 
			}
			if(*verbose)printf("lambda=%1.2E\tmaxdiff=%1.15f\n",lambda[ilambda],maxdiff);
			R_CheckUserInterrupt();
		}while(maxdiff*nzero>*eps);
	}
	free(resid);
}


#define oneThird 0.3333333333333333

void p7(double *x, double *y, int *n, int *p, double *beta, 
				double *lambda, double *alpha, int *nlambda, double *eps, int *niter, int *verbose)
// pearson type vii prior (scaled t)
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ilambda, thislambdaN, iter;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, aa, sqrtK, theta, root1, root2, root3;

	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=0;
	for(i=0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
//		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];


	for(ilambda=(*nlambda)-1; ilambda>=0; --ilambda) {
		thislambdaN=ilambda*(*p);
		aa=alpha[ilambda]*alpha[ilambda];
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nlambda)), sizeof(double)*(*p)); // hotstart
//printf("thislambdaN=%d\tlastlambdaN=%d\n",thislambdaN,  ((thislambdaN+*p)%(*p**nlambda)));continue;
		iter=0;
		do{
			maxdiff=negLarge; //nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
                A=-lsb; B=2.0*lambda[ilambda]+aa; C=-lsb*aa;
                Mterm=2.0*A*A*A-9.0*A*B+27.0*C; K=A*A-3.0*B; Nterm=(Mterm*Mterm-4.0*K*K*K);
				lastbeta=beta[thislambdaN+j];
/*  There is no need to check convexity, as long as global minimum is used
#define convexCheck(x) ((((x)*(x)+aa-lambda[ilambda])*((x)*(x)+aa-lambda[ilambda])-lambda[ilambda]*lambda[ilambda]+4*lambda[ilambda]*aa)>=0? x :beta[thislambdaN+j])
*/
				if(Nterm>0.0){
					sqrtNterm=sqrt(Nterm);

					beta[thislambdaN+j]=
						-oneThird*(A+R_pow(fabs(Mterm+sqrtNterm)*.5, oneThird)*sign(Mterm+sqrtNterm)+
									 R_pow(fabs(Mterm-sqrtNterm)*.5, oneThird)*sign(Mterm-sqrtNterm));
				}else{  // three real roots
					
					sqrtK=sqrt(K);
					theta=acos(-.5*Mterm/sqrtK/K);
					root1=2.0*oneThird*sqrtK*cos(theta*oneThird)+lsb*oneThird;
					root2=2.0*oneThird*sqrtK*cos((M_2PI+theta)*oneThird)+lsb*oneThird;
					root3=2.0*oneThird*sqrtK*cos((2.0*M_2PI+theta)*oneThird)+lsb*oneThird;

					root1=.5*(lsb-root1)*(lsb-root1)+lambda[ilambda]*log1p(root1*root1/aa) <
						  .5*(lsb-root2)*(lsb-root2)+lambda[ilambda]*log1p(root2*root2/aa) ?
						  root1 : root2 ;
					beta[thislambdaN+j]=
						 .5*(lsb-root1)*(lsb-root1)+lambda[ilambda]*log1p(root1*root1/aa) <
						 .5*(lsb-root3)*(lsb-root3)+lambda[ilambda]*log1p(root3*root3/aa) ?
						 root1 : root3 ;
				}
//printf("***\tj=%d\tlast.j=%d\tlastjN=\%d\tthisjN=%d\n",j,lastj, lastjN,thisjN);
//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
//for(i=0;i<*n;++i)printf("r[%d]=%f\t",i,resid[i]); printf("\n");
//				if(*verbose)printf("thisdiff=%1.15f\t",beta[thislambdaN+j]-lastbeta);
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 
			}
//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
			if(*verbose)printf("lambda=%1.2E\talpha=%1.2E\tmaxdiff=%1.15f\n",lambda[ilambda],alpha[ilambda],maxdiff);
			R_CheckUserInterrupt();
		}while(maxdiff>*eps && ++iter<*niter);
	}
	free(resid);
}

#define Factor (1.9802913004322129507587985355868732882822056101876)
#define log1pff (1.5936242600400400923230418758751602417890024248188)

void lassot(double *x, double *y, int *n, int *p, double *beta, 
				double *lambda, double *alpha, int *nlambda, double *eps, int *niter, int *verbose)
// pearson type vii prior for |x|<alpha (i.e. bayesA); lasso type prior for |x|>alpha
{	
	double thisdiff, maxdiff, maxreldiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ilambda, thislambdaN, lastlambda, iter;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, aa, lastaa, 
			sqrtK, theta, root1, root2, root3, tmproot,
			alphaFactor, lambdaprime;

	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=0;
	for(i=0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
//		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];


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




struct pMbeta
{	double * beta;
	double m;
	int p;
	double log1pbbOveraa;
};

double objForAlpha2(double alpha2, struct pMbeta * args)
{
	int j;
	double ans=0.0;
	for(j=0; j<(args->p); ++j)
		ans+=log1p(*(args->beta+j)**(args->beta+j)/alpha2);
	ans*=args->m;
	ans+=((double) args->p)*(lbeta(args->m-.5, .5));
	return ans;
}

void p7pmle(double *x, double *y, int *n, int *p, double *beta, 
				double *alpha, double *eps, int *verbose)
// pearson type vii prior; with all sigma^2e, alpha, and M updated in each cycle.
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ilambda, thislambdaN;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, sqrtK, theta, root1, root2, root3;
	double RSS, bbOveraabb, aa, sig2e, lambda, lastaa, lastM, M; 
	struct pMbeta args;
	

	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=0;
	for(i=0,RSS=0.0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];

	sig2e=RSS/(*n);
	aa=*alpha**alpha;
	for(j=0, bbOveraabb=0.0; j<*p; ++j) bbOveraabb+=beta[j]*beta[j]/(aa+beta[j]*beta[j]);
	M=0.5*(*p)/bbOveraabb;
	lambda=sig2e*(M);
	
	args.beta=beta;
	args.p=*p;
	args.m=M;
	// iteration

	do{
		maxdiff=negLarge; //nzero=0;
		for(j=-1; j<(*p)-1; ){
			lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
			thisjN=j*(*n);
			for(i=0, lsb=0.0; i<*n; ++i){
				resid[i]+=*(x+thisjN+i)*beta[j]-*(x+lastjN+i)*beta[lastj];
				lsb+=*(x+thisjN+i)*resid[i];
			}
			A=-lsb; B=2.0*lambda+aa; C=-lsb*aa;
			Mterm=2.0*A*A*A-9.0*A*B+27.0*C; K=A*A-3.0*B; Nterm=(Mterm*Mterm-4.0*K*K*K);
			lastbeta=beta[j];
			if(Nterm>0.0){
				sqrtNterm=sqrt(Nterm);
				beta[j]=
					-oneThird*(A+R_pow(fabs(Mterm+sqrtNterm)*.5, oneThird)*sign(Mterm+sqrtNterm)+
								 R_pow(fabs(Mterm-sqrtNterm)*.5, oneThird)*sign(Mterm-sqrtNterm));
			}else{  // three real roots
				sqrtK=sqrt(K);
				theta=acos(-.5*Mterm/sqrtK/K);
				root1=2.0*oneThird*sqrtK*cos(theta*oneThird)+lsb*oneThird;
				root2=2.0*oneThird*sqrtK*cos((M_2PI+theta)*oneThird)+lsb*oneThird;
				root1=.5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
					  .5*(lsb-root2)*(lsb-root2)+lambda*log1p(root2*root2/aa) ?
					  root1 : root2 ;
				root3=2.0*oneThird*sqrtK*cos((2.0*M_2PI+theta)*oneThird)+lsb*oneThird;
				beta[j]=
					 .5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
					 .5*(lsb-root3)*(lsb-root3)+lambda*log1p(root3*root3/aa) ?
					 root1 : root3 ;
			}
			if(fabs(beta[j]-lastbeta)>maxdiff) 
				maxdiff=fabs(beta[j]-lastbeta) ; 

			RSS+=2.0*lsb*(lastbeta-beta[j])+beta[j]*beta[j]-lastbeta*lastbeta;
			bbOveraabb+=beta[j]*beta[j]/(aa+beta[j]*beta[j])-lastbeta*lastbeta/(aa+lastbeta*lastbeta);
//			printf("RSS=%f\tadj=%f\n",RSS,2.0*lsb*(lastbeta-beta[j])+beta[j]*beta[j]-lastbeta*lastbeta);
		}
//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
		if(*verbose)printf("lambda=%1.2E\tm=%1.2E\talpha^2=%1.2E\tmaxdiff=%1.15f\n",lambda,M, aa,maxdiff);
		R_CheckUserInterrupt();
		
		// update other parameters than betas for next iteration
		sig2e=RSS/(*n);
		M=0.5*(*p)/bbOveraabb;
		lambda=sig2e*(M);
		args.m=M;
		lastaa=aa;
		aa=Brent_fmin(1e-6, 1e6, (double (*)(double, void*)) objForAlpha2, &args, *eps);
		if(fabs(lastaa-aa)>maxdiff) maxdiff=fabs(lastaa-aa);
	
	}while(maxdiff>*eps);
	*alpha=sqrt(aa);
	free(resid);
}


double objForM(double M, struct pMbeta * args)
{
	double ans;
	ans=M*(args->log1pbbOveraa)+((double) (args->p))*lbeta(M-.5, .5);
	return ans;
}
void p7alpha(double *x, double *y, int *n, int *p, double *beta, 
				double *alpha, int *nalpha, double *eps, int *verbose)
// pearson type vii prior; alpha^2 treated as the tuning parameter
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ialpha, thislambdaN;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, sqrtK, theta, root1, root2, root3;
	double RSS, log1pbbOveraa, aa, sig2e, lambda, lastaa, lastM, M; 
	struct pMbeta args;
	

	// initialization
	resid   =(double *)malloc(sizeof(double)*(*n));
	thislambdaN=((*nalpha)-1)*(*p);
	for(i=0,RSS=0.0;i<*n;i++){
		resid[i]=y[i];
		for(j=0;j<*p;j++)resid[i]-=*(x+j*(*n)+i)*beta[thislambdaN+j];
		RSS+=resid[i]*resid[i];
	}
	for(i=0;i<*n;i++)resid[i]+=*(x+((*p)-1)*(*n)+i)*beta[thislambdaN+(*p)-1];

	sig2e=RSS/(*n);
	args.beta=beta;
	args.p=*p;
	
	// iteration

	for(ialpha=(*nalpha)-1; ialpha>=0; --ialpha) {
		thislambdaN=ialpha*(*p);
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nalpha)), sizeof(double)*(*p)); // hotstart
		aa=alpha[ialpha]*alpha[ialpha];
		for(j=0, log1pbbOveraa=0.0; j<*p; ++j) 
			log1pbbOveraa+=log1p(beta[thislambdaN+j]*beta[thislambdaN+j]/aa);
		args.log1pbbOveraa=log1pbbOveraa;
		M=Brent_fmin(.5+1e-6, 1e6, (double (*)(double, void*)) objForM, &args, *eps);
		lambda=sig2e*(M);
printf("\nnew.a^2=%1.2e\tnew.lambda=%1.2E\tnew.M=%1.2E\tRSS=%f\tlog1p(bb/aa)=%1.2E\n", aa,lambda, M, RSS,log1pbbOveraa);		
		do{
			maxdiff=negLarge; //nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
				A=-lsb; B=2.0*lambda+aa; C=-lsb*aa;
				Mterm=2.0*A*A*A-9.0*A*B+27.0*C; K=A*A-3.0*B; Nterm=(Mterm*Mterm-4.0*K*K*K);
				lastbeta=beta[thislambdaN+j];
				if(Nterm>0.0){
					sqrtNterm=sqrt(Nterm);
					beta[thislambdaN+j]=
						-oneThird*(A+R_pow(fabs(Mterm+sqrtNterm)*.5, oneThird)*sign(Mterm+sqrtNterm)+
									 R_pow(fabs(Mterm-sqrtNterm)*.5, oneThird)*sign(Mterm-sqrtNterm));
				}else{  // three real roots
					sqrtK=sqrt(K);
					theta=acos(-.5*Mterm/sqrtK/K);
					root1=2.0*oneThird*sqrtK*cos(theta*oneThird)+lsb*oneThird;
					root2=2.0*oneThird*sqrtK*cos((M_2PI+theta)*oneThird)+lsb*oneThird;
					root1=.5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
						  .5*(lsb-root2)*(lsb-root2)+lambda*log1p(root2*root2/aa) ?
						  root1 : root2 ;
					root3=2.0*oneThird*sqrtK*cos((2.0*M_2PI+theta)*oneThird)+lsb*oneThird;
					beta[thislambdaN+j]=
						 .5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
						 .5*(lsb-root3)*(lsb-root3)+lambda*log1p(root3*root3/aa) ?
						 root1 : root3 ;
				}
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 

				RSS+=(2.0*lsb*(lastbeta-beta[thislambdaN+j])+
						beta[thislambdaN+j]*beta[thislambdaN+j]-lastbeta*lastbeta);
				log1pbbOveraa+=log1p(beta[thislambdaN+j]*beta[thislambdaN+j]/aa)
							  -log1p(lastbeta*lastbeta/aa);
				
//				if(*verbose){
//					for(tmpi=0,rsstmp=0.0;tmpi<*n;tmpi++){
//						tmpy=0.0;
//						for(tmpj=0;tmpj<*p;tmpj++)tmpy+=*(x+tmpj*(*n)+tmpi)*beta[thislambdaN+tmpj];
//						rsstmp+=(y[tmpi]-tmpy)*(y[tmpi]-tmpy);
//					}
//					printf("rss=%f\tRSS=%f\tlsb=%f\tlast.beta=%f\tthis.beta=%f\tadj=%f\n",rsstmp,
//						RSS,lsb,lastbeta,beta[thislambdaN+j],
//						2.0*lsb*(lastbeta-beta[thislambdaN+j])+
//							beta[thislambdaN+j]*beta[thislambdaN+j]-lastbeta*lastbeta);
//				}
			}
	//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
			if(*verbose)printf("lambda=%1.2E\tm=%1.2E\talpha^2=%1.2E\tRSS=%f\tlog1p(bb/aa)=%1.2E\tmaxdiff=%1.15f\n",lambda,M, aa, RSS, log1pbbOveraa,maxdiff);
			R_CheckUserInterrupt();
			
			// update lambda, M and sig2e
			sig2e=RSS/(*n);
			args.log1pbbOveraa=log1pbbOveraa;
			M=Brent_fmin(.5+1e-6, 1e6, (double (*)(double, void*)) objForM, &args, *eps);
			lambda=sig2e*(M);
		
		}while(maxdiff>*eps);  // end for this alpha
	};  // end for the alpha loop
	free(resid);
}



void p7M(double *x, double *y, int *n, int *p, double *beta, 
				double *alpha, int *nalpha, double *eps, int *verbose)
// pearson type vii prior; M treated as the tuning parameter
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ialpha, thislambdaN;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, sqrtK, theta, root1, root2, root3;
	double RSS, log1pbbOveraa, aa, sig2e, lambda, lastaa, lastM, M; 
	struct pMbeta args;
	
	resid   =(double *)malloc(sizeof(double)*(*n));
	memcpy(resid,y,sizeof(double)*(*n));

	// initialization
	for(i=0, RSS=0.0; i<*n; ++i) RSS+=y[i]*y[i];
	sig2e=RSS/(*n);
	args.beta=beta;
	args.p=*p;
	
	// iteration

	for(ialpha=(*nalpha)-1; ialpha>=0; --ialpha) {
		thislambdaN=ialpha*(*p);
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nalpha)), sizeof(double)*(*p)); // hotstart
		aa=alpha[ialpha]*alpha[ialpha];
		for(j=0, log1pbbOveraa=0.0; j<*p; ++j) 
			log1pbbOveraa+=log1p(beta[thislambdaN+j]*beta[thislambdaN+j]/aa);
		args.log1pbbOveraa=log1pbbOveraa;
		M=Brent_fmin(.5+1e-6, 1e6, (double (*)(double, void*)) objForM, &args, *eps);
		lambda=sig2e*(M);
printf("\nnew.a^2=%1.2e\tnew.lambda=%1.2E\tnew.M=%1.2E\tRSS=%1.2E\tlog1p(bb/aa)=%1.2E\n", aa,lambda, M, RSS,log1pbbOveraa);		
		do{
			maxdiff=negLarge; //nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}
				A=-lsb; B=2.0*lambda+aa; C=-lsb*aa;
				Mterm=2.0*A*A*A-9.0*A*B+27.0*C; K=A*A-3.0*B; Nterm=(Mterm*Mterm-4.0*K*K*K);
				lastbeta=beta[thislambdaN+j];
				if(Nterm>0.0){
					sqrtNterm=sqrt(Nterm);
					beta[thislambdaN+j]=
						-oneThird*(A+R_pow(fabs(Mterm+sqrtNterm)*.5, oneThird)*sign(Mterm+sqrtNterm)+
									 R_pow(fabs(Mterm-sqrtNterm)*.5, oneThird)*sign(Mterm-sqrtNterm));
				}else{  // three real roots
					sqrtK=sqrt(K);
					theta=acos(-.5*Mterm/sqrtK/K);
					root1=2.0*oneThird*sqrtK*cos(theta*oneThird)+lsb*oneThird;
					root2=2.0*oneThird*sqrtK*cos((M_2PI+theta)*oneThird)+lsb*oneThird;
					root1=.5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
						  .5*(lsb-root2)*(lsb-root2)+lambda*log1p(root2*root2/aa) ?
						  root1 : root2 ;
					root3=2.0*oneThird*sqrtK*cos((2.0*M_2PI+theta)*oneThird)+lsb*oneThird;
					beta[thislambdaN+j]=
						 .5*(lsb-root1)*(lsb-root1)+lambda*log1p(root1*root1/aa) <
						 .5*(lsb-root3)*(lsb-root3)+lambda*log1p(root3*root3/aa) ?
						 root1 : root3 ;
				}
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 

				RSS+=2.0*lsb*(lastbeta-beta[thislambdaN+j])+
						beta[thislambdaN+j]*beta[thislambdaN+j]-lastbeta*lastbeta;
				log1pbbOveraa+=log1p(beta[thislambdaN+j]*beta[thislambdaN+j]/aa)
							  -log1p(lastbeta*lastbeta/aa);
			}
	//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
			if(*verbose)printf("lambda=%1.2E\tm=%1.2E\talpha^2=%1.2E\tRSS=%1.2E\tlog1p(bb/aa)=%1.2E\tmaxdiff=%1.15f\n",lambda,M, aa, RSS, log1pbbOveraa,maxdiff);
			R_CheckUserInterrupt();
			
			// update lambda, M and sig2e
			sig2e=RSS/(*n);
			args.log1pbbOveraa=log1pbbOveraa;
			M=Brent_fmin(.5+1e-6, 1e6, (double (*)(double, void*)) objForM, &args, *eps);
			lambda=sig2e*(M);
		
		}while(maxdiff>*eps);  // end for this alpha
	};  // end for the alpha loop
	free(resid);
}


void BC(double *x, double *y, int *n, int *p, double *pi0, double *beta, 
				double *alpha, int *nalpha, double *eps, int *verbose)
// BayesC prior; update to the posterior mean in each iterations; alpha is the ratio of sig2e/sig2b
{	
	double thisdiff, maxdiff;
	double lastbeta, *resid;
	int i,j, thisjN, lastj, lastjN, ialpha, thislambdaN;
	double lsb, A, B, C, K, Mterm, Nterm, sqrtNterm, sqrtK, theta, root1, root2, root3;
	double RSS, log1pbbOveraa, aa, sig2e, lambda, lastaa, lastM, M; 
	struct pMbeta args;
	
	resid   =(double *)malloc(sizeof(double)*(*n));
	memcpy(resid,y,sizeof(double)*(*n));

	// initialization
	for(i=0, RSS=0.0; i<*n; ++i) RSS+=y[i]*y[i];	
	if(*verbose)printf("RSS=%1.2E\n", RSS);
	// iteration

	for(ialpha=(*nalpha)-1; ialpha>=0; --ialpha) {
		thislambdaN=ialpha*(*p);
		memcpy(beta+thislambdaN, beta+((thislambdaN+*p)%(*p**nalpha)), sizeof(double)*(*p)); // hotstart
		do{
			maxdiff=negLarge; //nzero=0;
			for(j=-1; j<(*p)-1; ){
				lastjN=(lastj=(*p+(j++))%(*p))*(*n); // modulo is not safe for negative numbers
				thisjN=j*(*n);
				for(i=0, lsb=0.0; i<*n; ++i){
					resid[i]+=*(x+thisjN+i)*beta[thislambdaN+j]-*(x+lastjN+i)*beta[thislambdaN+lastj];
					lsb+=*(x+thisjN+i)*resid[i];
				}

				lastbeta=beta[thislambdaN+j];
				beta[thislambdaN+j]=(1.0-*pi0)*lsb/(1.0+alpha[ialpha]);
				if(fabs(beta[thislambdaN+j]-lastbeta)>maxdiff) 
					maxdiff=fabs(beta[thislambdaN+j]-lastbeta) ; 

				RSS+=(2.0*lsb*(lastbeta-beta[thislambdaN+j])+
						beta[thislambdaN+j]*beta[thislambdaN+j]-lastbeta*lastbeta);
			}
	//printf("lsb=%f\tA=%f\tB=%f\tC=%f\tM=%f\tK=%f\tN=%f\n",lsb,A,B,C,Mterm,K,sqrtNterm); //return ;
			if(*verbose)printf("alpha=%1.2E\tRSS=%1.2E\tmaxdiff=%1.15f\n",alpha[ialpha], RSS,maxdiff);
			R_CheckUserInterrupt();
			
		}while(maxdiff>*eps);  // end for this alpha
	};  // end for the alpha loop
	free(resid);
}
