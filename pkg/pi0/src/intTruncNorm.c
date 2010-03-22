#ifdef _WINDOWS
#define INLINE __inline
#ifdef INTNCT_EXPORTS
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif
#else 
#define DLL_API
#define INLINE inline
#endif
 
#include<Rmath.h>
#include<R_ext/Arith.h>
#include<stdlib.h>
#ifdef _DEBUG
#include<stdio.h>
#endif

DLL_API INLINE void intTruncNorm(int *r, int *length, double *mu, double *sd, double *low, double *upp, double *ans)
{
	// r needs to be non-negative; length needs to be positive
	// returns int_low^upp x^(r-1+1:length) exp(-(x-mu)^2/2/sd^2) dx

	double sd2=*sd*(*sd);

	double 	expUpp=sd2*exp(-(*upp-*mu)*(*upp-*mu)*.5/sd2);   // NAOK
	double 	expLow=sd2*exp(-(*low-*mu)*(*low-*mu)*.5/sd2);	// NAOK

	double 	Ir_2=2.506628274631*(*sd)*(pnorm(*upp, *mu, *sd, 1, 0)-pnorm(*low, *mu, *sd, 1, 0));	// NAOK
	double  Ir_1=(*mu)*Ir_2 + (expLow-expUpp);
#ifdef _DEBUG
	printf("\n\nsd^2=%f\texpUpp=%f\texpLow=%f\tIr_2=%f\tIr_1=%f\tpnormUpp=%f\tpnormLow=%f\n",
		sd2,expUpp,expLow,Ir_2,Ir_1, pnorm(*upp, *mu, *sd, 1, 0), pnorm(*low, *mu, *sd, 1, 0));
#endif
	double Ir;
	double lastUppPow=(*upp);
	double lastLowPow=(*low);

	int finiteUpp=R_finite(lastUppPow);
	int finiteLow=R_finite(lastLowPow);
#ifdef _DEBUG
	printf("finiteUpp=%d\tfiniteLow=%d\n",finiteUpp, finiteLow) ;
#endif

	int out=0;
	int i;

	if((*r)==0) ans[out++]=Ir_2;
	if((*r)==1 || (*length>1 && *r==0) ) ans[out++]=Ir_1;

	for(i=2; i<=(*r)+(*length)-1; ++i){
#ifdef _DEBUG
	printf("i=%d\tIr_1=%f\tIr_2=%f\tlastUppPow finite=%d\tlastLowPow finite=%d\tadditive=%f\n", 
		i, Ir_1, Ir_2, R_finite(lastUppPow), R_finite(lastLowPow)
		, - (finiteUpp ? lastUppPow*expUpp : 0.0) + (finiteLow ? lastLowPow*expLow : 0.0)
	) ;
#endif
		Ir = sd2*(i-1)*Ir_2 + (*mu)*Ir_1 - (finiteUpp ? lastUppPow*expUpp : 0.0) + (finiteLow ? lastLowPow*expLow : 0.0) ; // deal NA
		Ir_2=Ir_1;
		Ir_1=Ir;

		lastUppPow*=(*upp);	// NAOK
		lastLowPow*=(*low);	// NAOK
		if(i>=*r) ans[out++]=Ir;
	}
}


DLL_API INLINE double divDifInt(int ntab, int r, double *ytab, double xval, int takeLog)	// warning: this function will modify values of ytab
{	// reference: http://people.sc.fsu.edu/~burkardt/cpp_src/divdif/divdif.html
	int i,j;
	double value;

//	for(i=0;i<ntab;++i) printf("%f\t", ytab[i]); printf("\n");
	if(takeLog)	for(i=0;i<ntab;++i) if(ytab[i]<=0) {takeLog =0; break;}
	if(takeLog)	for(i=0;i<ntab;++i) ytab[i]=log(ytab[i]);
#ifdef _DEBUG
	for(i=0; i<ntab; ++i) printf("ytab[%d]=%f\t", i, ytab[i]); printf("\n");
#endif
	
	for(i=1; i<= ntab-1; ++i)
		for(j=ntab-1; i<=j; --j)
			ytab[j]=(ytab[j]-ytab[j-1])/ i ;

	value=ytab[ntab-1];
#ifdef _DEBUG
	printf("value(before eval)=%f\n", value);
#endif

	for(i=2; i<=ntab; ++i)
		value=ytab[ntab-i] + (xval-(r+ntab-i) )*value;
#ifdef _DEBUG
	printf("value(after eval)=%f\n", value);
#endif

	return takeLog? exp(value) : value; 
}

DLL_API INLINE void fracTruncNorm(double *r, double *mu, double *sd, double *low, double *upp, double *ans, int *ndiv, int* takeLog)
{
	// r needs to be non-negative; ndiv needs to be positive
	// returns int_low^upp x^r exp(-(x-mu)^2/2/sd^2) dx

	double *neighbors=(double*)malloc(*ndiv*sizeof(double));
	int integerR=ceil(*r-*ndiv*.5) >0 ? (int)ceil(*r-*ndiv*.5) : (int)0;

	intTruncNorm(&integerR, ndiv, mu, sd, low, upp, neighbors);
#ifdef _DEBUG
	printf("ndiv=%d\tintegerR=%d\tr=%f\nneighbors[0]=", *ndiv, integerR, *r);
	int i;
	for(i=0;i<*ndiv;++i) printf("%f\t",neighbors[i]); printf("\n");
#endif
	*ans=divDifInt(*ndiv, integerR, neighbors, *r, *takeLog);
	free(neighbors);
}


DLL_API void intTruncNormVec(int *n, int *r, double *mu, double *sd, double *low, double *upp, double *ans)
{	// all vectors: r, mu, sd, low, upp, ans should be of length *n
	int i; 
	int len=1;
#ifdef _DEBUG
	for(i=0; i<*n; ++i) printf("n=%d\tr=%d\tmu=%f\tsd=%f\t%low=%f\tupp=%f\tans=%f\n",
							  *n, *(r+i), *(mu+i), *(sd+i), *(low+i), *(upp+i), *(ans+i));
#endif
	for(i=0; i<*n; ++i) intTruncNorm(r+i, &len, mu+i, sd+i, low+i, upp+i, ans+i);
}

DLL_API void fracTruncNormVec(int *n, double *r, double *mu, double *sd, double *low, double *upp, double *ans, int *ndiv, int* takeLog)
{	// all vectors: r, mu, sd, low, upp, ans should be of length *n
	int i;
#ifdef _DEBUG
	for(i=0; i<*n; ++i) printf("n=%d\tr=%f\tmu=%f\tsd=%f\tlow=%f\tupp=%f\tans=%f\tndiv=%d\ttakeLog=%d\n",
							  *n, *(r+i), *(mu+i), *(sd+i), *(low+i), *(upp+i), *(ans+i), *ndiv, *takeLog);
#endif
	for(i=0; i<*n; ++i) fracTruncNorm(r+i, mu+i, sd+i, low+i, upp+i, ans+i, ndiv, takeLog);
}
