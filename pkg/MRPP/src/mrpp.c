#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#undef DEBUG

#ifdef DEBUG
	int debug_i, debug_j, debug_k, debug_n, debug_m, debug_N;
	
	void printIntVec(int *x, int n, char * head)
	{
		Rprintf("%s ", head);
		for(debug_i=0; debug_i<n; debug_i++)
			Rprintf("%3d", x[debug_i]);
		Rprintf("\n");
	}
#endif

extern void radixsort(int *, int , int *);
/*
//void unisqeuc(double *x, int *nr,  double *dist)
*********** CALL FROM R ***************
.C("unisqeuc", x = as.double(1:2326), nr = as.integer(2326), dist=dist.holder, DUP = FALSE,  NAOK = TRUE)$dist 
    where x is univariate double vector; 
         nr is a scalar integer equal to the length x; 
         dist is a double vector of length nr*nr, which would store the resulting squared euclidean distance matrix
***************************************
//{
//    int i, j, basei, basej, dimi;
//    for(j = 1, basej=*nr ; j < *nr ; ++j, basej+=*nr){
//            for(i = basei=0 ; i < j ; ++i, basei+=*nr)
//                dist[i*(*nr)+j]=dist[basej+i] += (*(x+j)-*(x+i))*(*(x+j)-*(x+i));
//    }
//}


//void sqeuc(double *x, int *nr, int *nc, double *dist)
//{
//    double *x0;
//    int i, j, basei, basej, dimi;
//    for(dimi=0; dimi<*nc; ++dimi){
//        x0=x+dimi**nr;
//        for(j = 1, basej=*nr ; j < *nr ; ++j, basej+=*nr){
//            for(i = basei=0 ; i < j ; ++i, basei+=*nr)
//                dist[basej+i] += (*(x0+j)-*(x0+i))*(*(x0+j)-*(x0+i));
//        }
//    }
//    for(j=1, basej=*nr; j<*nr;++j, basej+=*nr)
//        for(i=basei=0;i<j;++i,basei+=*nr)
//            dist[i*(*nr)+j]=dist[basej+i];
//}

    
//#include<stdio.h>

//inline double sumSubMat( double *x,  int *idx,  int n,  int N)
//// Input:
////   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
////     idx is a vector of indices where sum of x is taken over;
////   n is the length of idx; 
////   *N is the total sample size; 
//// Output: 
////   returns sum(as.matrix(x)[idx, idx])
//{
//    double ans=0.0;
//    register unsigned  int i,j,row,col, N2=(N<<1);
//
//    for(i=n-1; i; --i){
//        for(j=0; j<i; ++j){
////            if(idx[i]>idx[j]) row=idx[i],col=idx[j]; 
////            else (row=idx[j],col=idx[i]);
////            ans+=x[ ((col-1)*N-((col*col-col)>>1)+row-col) -1];
//                // column head index (R) is ((col-1)*N-((col*col-col)>>1)+1
//                // row adjustment is  row-col-1
//                // C indexing adjustment is -1
//                // Note that >> has lower precedence than +/-
//            idx[i]>idx[j] ? (row=idx[i],col=idx[j]) : (row=idx[j],col=idx[i]);
//            ans+=x[ ((col*(N2-col-1))>>1)-N+row-1 ]; // equivalent! but faster?
//        }
//    }
//    ans*=2;
//    return ans;
//}
*/

static R_INLINE double sumSubMatSorted( double const * const x,  int const * const idx,  const int n,  const int N)
/*
// Input:
//   x is a stacked vector of the lower triangle of distance mat of (N x N); 
//     idx is a vector of "sorted" indices over which sum of x is taken (idx starts from 1 instead of 0);
//   n is the length of idx; 
//   N is the total sample size; 
// Output: 
//   returns sum(as.matrix(x)[idx, idx])
// Note:
//   This is about 3~5% faster than sumSubMat when N=30 and B=5000. 
*/
{
    double ans, compen, adjx, tmpAns;
    register unsigned int j, i, N2=(N<<1);
    register int base;

    ans = compen = 0.0;
    for(j=0; j<n-1; ++j){ /* // column index */
        /*  base=((idx[j]*(N2-idx[j]-1))>>1)-N-1;   // part that does not involve row index */
        base = (((N2-idx[j])*(idx[j]-1))>>1) - idx[j] - 1 ;  /* this should replace the previous line  */
        for(i=j+1; i<n; ++i){  /* //  row index */
            /*  ans+= x [base + idx[i]];   // the following Kahan's algo recovers this line */
            adjx = x[base + idx[i]] - compen;
            tmpAns = ans + adjx;
            compen = (tmpAns - ans) - adjx;
            ans = tmpAns;
        }
    }
/*    ans*=2.0;    */
    return ans*2.0;
}

#ifdef DEBUG
void testSumSubMatSorted(double const * const x, int const * const idx, const int * const n, const int * const N, double  * const ans)
{
    *ans=sumSubMatSorted(x, idx, *n, *N);
}

void mrppstats2(double const * const x, int const * const perm, int const * const permcomplement, 
               int const * const n, int const * const B, int const * const N, int const * const wtmethod,
               double * const ans)
{
/*
// Input:
//   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
//   perm is perm idx mat for 1st group (*n x *B); 
//   permcomplement is perm idx mat for 2nd group ((*N-*n) x *B); 
//   *n is the number of rows in perm; *B is the number of cols in perm; 
//   *N is the total sample size; 
//   *wtmethod is the treatment group weight: 0=sample size-1; 1=sample size
// Output: 
//   *ans is a length *B space of mrpp statistics to be filled in
*/
    register unsigned int b, Nn;
    double tmp1, tmp2, denom1, denom2;

    Nn = *N - *n;
    denom1=1.0/(*n-*wtmethod); denom2=1.0/(Nn-*wtmethod) ;
    for(b=0; b<*B; ++b) 
    {
/*
//        tmp1=sumSubMat(x, perm+(b**n), *n, *N); 
//        tmp2=sumSubMat(x, permcomplement+(b*(*N-*n)), *N-*n, *N); 
*/
        tmp1=sumSubMatSorted(x, perm+(b**n), *n, *N); 
        tmp2=sumSubMatSorted(x, permcomplement+(b*Nn), Nn, *N); 
        

           ans[b]=tmp1*denom1 +  tmp2*denom2 ;
    }
}

#endif

	
static R_INLINE SEXP mrppstats_listOfMatrix(double * ptrY, SEXP permMats, double wt, R_len_t N)
{
	/* ptrY is the REAL pointer to the vector of a 'dist' object with each element being double */
	/* permMats is a 'list', of n * B permutation indices; see permuteTrt. */
	/* wt is a double scaler of 0 or 1 weighting method: 0.0=sample size-1; 1.0=sample size */
	/* N is the total sample size */
    SEXP ans;
    R_len_t t, B, ntrt, b, n;
    double  * ptrAns, denom, dn;
	int * ptrPerm;
    
    ntrt=LENGTH(permMats);
    B=Rf_ncols(VECTOR_ELT(permMats, 1));

    PROTECT(ans = NEW_NUMERIC(B));
    ptrAns=REAL(ans);
	for(b=0; b<B; ++b) *(ptrAns++)=0.0;
	ptrAns=REAL(ans);

#ifdef DEBUG
	Rprintf("ntrt=%d\tB=%d\tN=%d\twt=%f\n", ntrt, B, N, wt);
#endif

    for(t=ntrt-1; t>=0; --t){
        n = Rf_nrows(VECTOR_ELT(permMats, t));
        dn = (double) n;
        denom = 1.0 / (dn - wt) ;
#ifdef DEBUG
		Rprintf("t=%d\tn=%d\tdn=%f\tdenom=%f\n", t, n, dn, denom);
#endif
        ptrPerm = INTEGER(VECTOR_ELT(permMats, t)); 
        for(b=0; b<B; ++b){
            ptrAns[b] += sumSubMatSorted(ptrY,  ptrPerm + (b * n)  , n, N) * denom;
        }
    }
    UNPROTECT(1);
    return(ans);
}

static R_INLINE SEXP mrppstats_string(double * ptrY, SEXP permString, SEXP perm0, double wt, R_len_t N)
{
	/* ptrY is the REAL pointer to the vector of a 'dist' object with each element being double */
	/* permMats is a permutedTrt object with 'idx' attribute being a character vector. */
	/* wt is a double scaler of 0 or 1 weighting method: 0.0=sample size-1; 1.0=sample size */
	/* N is the total sample size */
    SEXP ans, perm, NSEXP, string1;
    R_len_t j, t, B, ntrt, b, *ns, **ptrPerm0, *permBuf, *permBuf0, *sortBuf, * ptrPerm;
    double  * ptrAns,  dn, *denoms;
	SEXP dec2pvCall, dec2pvCall2;
	 
	PROTECT(NSEXP = NEW_INTEGER(1));
	*(INTEGER(NSEXP)) = N;
	
    ntrt=LENGTH(perm0);  // LENGTH(permMats);
    B=LENGTH(permString) ;  // Rf_ncols(VECTOR_ELT(permMats, 1));

    PROTECT(ans = NEW_NUMERIC(B));
    ptrAns=REAL(ans);
	for(b=0; b<B; ++b) *(ptrAns++)=0.0;
	ptrAns=REAL(ans);

#ifdef DEBUG
	Rprintf("ntrt=%d\tB=%d\tN=%d\twt=%f\n", ntrt, B,  N, wt);
#endif

	ns = (int *) R_alloc(ntrt, sizeof(int)); 
	denoms = (double *) R_alloc(ntrt, sizeof(double)); 
	ptrPerm0 = (int **) R_alloc(ntrt, sizeof(int *)); 
	permBuf0 = permBuf = (int *) R_alloc(N, sizeof(int));
	sortBuf = (int *) R_alloc(N, sizeof(int));
	
	for(t=0; t < ntrt; ++t){
		ns[t] = LENGTH(VECTOR_ELT(perm0, t));
		dn = (double) ns[t];
        denoms[t] = 1.0 / (dn - wt) ;
		ptrPerm0[t] = INTEGER(VECTOR_ELT(perm0, t)); 
#ifdef DEBUG
		Rprintf("ns[%d]=%d\tptrPerm0[%d]=%d\n", t, ns[t], t, *(ptrPerm0[t]));
#endif
	}


	
	PROTECT(string1 = NEW_STRING(1));
	PROTECT(dec2pvCall = dec2pvCall2 = allocList(3));
	SET_TYPEOF(dec2pvCall, LANGSXP);
		SETCAR(dec2pvCall, install("dec2permvec"));  // this is a call. 
	dec2pvCall2 = CDR(dec2pvCall);
		SETCAR(dec2pvCall2, string1);				 // this string is not bound to any name but a literal
	dec2pvCall2 = CDR(dec2pvCall2);
		SETCAR(dec2pvCall2,  NSEXP); 				 // this is a literal too. 
#ifdef DEBUG
		SET_STRING_ELT(string1, 0, STRING_ELT(permString, 1));
		// UNPROTECT(4);
		// return dec2pvCall;
#endif

	for(b=0; b<B; ++b) {
		SET_STRING_ELT(string1, 0, STRING_ELT(permString, b));
		PROTECT(perm = EVAL(dec2pvCall));   // /* R:  perm=dec2permvec(decfrCC[b],N) */  // EVAL'ed in global env, because arguments are all literals. It does not matter where to eval this call. 
		ptrPerm = INTEGER(perm); 
#ifdef DEBUG
		Rprintf("b=%d\n\t", b);
		printIntVec(ptrPerm, N, "Permvec:");
#endif
		
		for(t=ntrt-1; t>=0; --t) {
#ifdef DEBUG
			Rprintf("\tt=%d\tn=%d\n\t",  t, ns[t] );
#endif
			// /* R: for(i in seq(ntrts)) ans[[i]][,b]=.Call(radixSort_prealloc, perm[part0[[i]]], buff) */
			for(j=0, permBuf=permBuf0; j<ns[t]; ++j) *(permBuf++) = *(ptrPerm + * (ptrPerm0[t]+j) - 1);
#ifdef DEBUG
			printIntVec(permBuf0, ns[t], "Pre-sort  :"); 
#endif	
			radixsort(permBuf0, ns[t] , sortBuf);
#ifdef DEBUG
			printIntVec(permBuf0, ns[t], "PermMatCol:"); 
#endif			
			ptrAns[b] += sumSubMatSorted(ptrY,  permBuf0, ns[t], N) * denoms[t];
		}
		
		UNPROTECT(1); // of perm
	 }

    UNPROTECT(4);
    return(ans);
}

SEXP mrppstats(SEXP y, SEXP Perms, SEXP wtmethod)
{
	/* y is the vector of a 'dist' object with each element being double */
	/* permMats is a 'permutedTrt' object. See permutedTrt R function */
	/* wtmethod is a scaler of 0 or 1 weighting method: 0=sample size-1; 1=sample size */
	R_len_t N;
	double  wt, *ptrY;
	SEXP wt_real, permString;
	
	if(isReal(wtmethod)) {
		wt = *(REAL(wtmethod));
	}else{
		PROTECT(wt_real = AS_NUMERIC(wtmethod));
		wt = *(REAL(wt_real));
		UNPROTECT(1);
	}

	ptrY = REAL(y);
    N = LENGTH(y);  // length of 'dist' obj
    N = ( 1 + (int)(0.5 + sqrt(1.0 + 8.0 * N)) ) >> 1;

	permString = getAttrib(Perms, install("idx")) ;
	if (R_NaString == STRING_ELT(permString, 0) ) {
		return mrppstats_listOfMatrix(ptrY, Perms, wt, N);
	}else {	
		return mrppstats_string(ptrY, permString, Perms, wt, N);
	}
}

SEXP sumThresh0(SEXP x)
// computing sum(pmax(x,0)) assuming x is a double vector;
// For speed, no type check is conducted! 
// This is no longer used in the R funciton smrpp.penWt. 
{
	SEXP ans;
	R_len_t   i;
	double * ptrAns, * ptrX; 
	 
	PROTECT( ans = NEW_NUMERIC(1) );
	ptrAns = REAL(ans);
	ptrX   = REAL(x);
	
	i=LENGTH(x); 
	*ptrAns = 0.0;
    for(; i>0; --i, ++ptrX) 
		if (*ptrX > 0.0) (*ptrAns) += *ptrX ;
	
	UNPROTECT(1);
	return ans;
}

SEXP objSolveDelta(SEXP del, SEXP dpdw, SEXP fact, SEXP b, SEXP l, SEXP R)
// computing fact[b,l]*sum(pmax(x,0)) - R;
// Equivalent R function is:    f=function(del) fact[b,l]*sum(pmax(-del-dp.dw[b, ], 0)) - R  
// b and l are integers; others are double;
// For speed, no type check is conducted! 
// This is only used in the old R funciton smrpp.penWt. (no longer used)
{
	double negDel;
	SEXP ans;
	R_len_t   i, B;
	double * ptrAns, * ptrX; 
	 
	PROTECT( ans = NEW_NUMERIC(1) );
	ptrAns = REAL(ans);
	ptrX   = REAL(dpdw) + *INTEGER(b) - 1;
	
	negDel = - (* REAL(del)); 

	i=Rf_ncols(dpdw); 
	B=Rf_nrows(dpdw);

	*ptrAns = 0.0;
    for(; i>0; --i, (ptrX += B) ) 
		if (*ptrX < negDel) (*ptrAns) += negDel - (*ptrX) ;
	
	*ptrAns = *ptrAns * (*(REAL(fact) + B * (*(INTEGER(l))-1) + *INTEGER(b) - 1)) - *REAL(R) ;
	UNPROTECT(1);
	return ans;
}


SEXP pmax0(SEXP x)
// computing (pmax(x,0)) assuming x is a double vector;
// For speed, no type check is conducted! 
// This is no longer used in the R funciton smrpp.penWt. 
{
	SEXP ans;
	R_len_t   i;
	double * ptrAns, * ptrX; 
	 
	i=LENGTH(x); 
	PROTECT( ans = NEW_NUMERIC(i) );
	ptrAns = REAL(ans);
	ptrX   = REAL(x);
	
    for(; i>0; --i, ++ptrX, ++ptrAns) 
		if (*ptrX > 0.0) (*ptrAns) = *ptrX ;
	    else (*ptrAns) = 0.0;
	
	UNPROTECT(1);
	return ans;
}
