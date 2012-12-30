#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP FR2permvec(SEXP FR, SEXP base)
{
	SEXP ans, cand/* , oFR, oBase */;
	R_len_t n, i, b, idx;
	int * intPtr, *headPtr, *FRPtr, *ansPtr; 

/* 	if (!IS_INTEGER(FR)) {
		oFR = FR;
		PROTECT(FR = AS_INTEGER(FR));
	}
	if (!IS_INTEGER(base)) {
		oBase = base;
		PROTECT(base = AS_INTEGER(base));
	} */
	
	b = *(INTEGER(base)); 
	n = LENGTH(FR);
	PROTECT(ans = NEW_INTEGER(n));
	
	PROTECT(cand = NEW_INTEGER(n));
	intPtr = headPtr = INTEGER(cand); 
	for(i=0; i<n; ++i) *(intPtr++) = i + b;
	
	FRPtr = INTEGER(FR); 
	ansPtr = INTEGER(ans);
	for(--n; n>=0; --n){
		idx = *(FRPtr++);
		*(ansPtr++) = *(headPtr + idx);
		// Rprintf("idx=%d; n=%d ; ", idx, n);
		if ((idx<<1) < n){
			memmove (headPtr + 1, headPtr, idx * sizeof(int) );
			++headPtr;
			// for(i=0; i<n; ++i) Rprintf("%d ", *(headPtr+i));
			// Rprintf("h\n");
		}else{
			memmove (headPtr + idx, headPtr + idx + 1, (n - idx) * sizeof(int));
			// for(i=0; i<n; ++i) Rprintf("%d ", *(headPtr+i));
			// Rprintf("t\n");
		}
	}
	UNPROTECT(2);
/* 	if (!IS_INTEGER(oBase)) UNPROTECT(1);	
	if (!IS_INTEGER(oFR)) UNPROTECT(1);	
 */	return(ans);
}
