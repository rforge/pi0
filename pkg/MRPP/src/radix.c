#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define RADIX (8) 

void radixsort(int *a, int n, int *b)
{
  static int bucket[RADIX]; 
  int i, m = a[0], exp = 1;
  for (i = 0; i < n; ++i)
  {
    if (a[i] > m)
      m = a[i];
  }
  
  
  while (m / exp > 0)
  {
	for(i=0; i<RADIX; ++i) bucket[i]=0;  // memset (bucket, 0, RADIX); 
    // int bucket[RADIX] =  {  0 };
    for (i = 0; i < n; ++i)
      bucket[a[i] / exp % RADIX]++;
    for (i = 1; i < RADIX; ++i)
      bucket[i] += bucket[i - 1];
    for (i = n - 1; i >= 0; --i)
      b[--bucket[a[i] / exp % RADIX]] = a[i];
    for (i = 0; i < n; i++)
      a[i] = b[i];
    exp *= RADIX;
   }
}

SEXP radixSort_prealloc(SEXP x, SEXP buff)
// Radix sort using pre-allocated buffer space. 
// For speed, this does not check any assertions.
{
	int * ptrAns;
	R_len_t n;
	SEXP ans;
	
	n=LENGTH(x);
	
	PROTECT(ans = NEW_INTEGER(n));
	ptrAns = INTEGER(ans);
	memcpy (ptrAns, INTEGER(x), sizeof(int) * n);
	
	radixsort(ptrAns, n, INTEGER(buff));
	
	UNPROTECT(1);
	return ans;
}