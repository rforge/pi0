/*
//void unisqeuc(double *x, int *nr,  double *dist)
*********** CALL FROM R ***************
.C("unisqeuc", x = as.double(1:2326), nr = as.integer(2326), dist=dist.holder, DUP = FALSE,  NAOK = TRUE)$dist 
	where x is univariate double vector; 
	     nr is a scalar integer equal to the length x; 
		 dist is a double vector of length nr*nr, which would store the resulting squared euclidean distance matrix
***************************************
//{
//	int i, j, basei, basej, dimi;
//	for(j = 1, basej=*nr ; j < *nr ; ++j, basej+=*nr){
//			for(i = basei=0 ; i < j ; ++i, basei+=*nr)
//				dist[i*(*nr)+j]=dist[basej+i] += (*(x+j)-*(x+i))*(*(x+j)-*(x+i));
//	}
//}


//void sqeuc(double *x, int *nr, int *nc, double *dist)
//{
//	double *x0;
//    int i, j, basei, basej, dimi;
//	for(dimi=0; dimi<*nc; ++dimi){
//		x0=x+dimi**nr;
//	    for(j = 1, basej=*nr ; j < *nr ; ++j, basej+=*nr){
//			for(i = basei=0 ; i < j ; ++i, basei+=*nr)
//				dist[basej+i] += (*(x0+j)-*(x0+i))*(*(x0+j)-*(x0+i));
//		}
//	}
//	for(j=1, basej=*nr; j<*nr;++j, basej+=*nr)
//		for(i=basei=0;i<j;++i,basei+=*nr)
//			dist[i*(*nr)+j]=dist[basej+i];
//}

	
//#include<stdio.h>

//inline double sumSubMat( double *x,  int *idx,  int n,  int N)
//// Input:
////   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
////	 idx is a vector of indices where sum of x is taken over;
////   n is the length of idx; 
////   *N is the total sample size; 
//// Output: 
////   returns sum(as.matrix(x)[idx, idx])
//{
//	double ans=0.0;
//	register unsigned  int i,j,row,col, N2=(N<<1);
//
//	for(i=n-1; i; --i){
//		for(j=0; j<i; ++j){
////			if(idx[i]>idx[j]) row=idx[i],col=idx[j]; 
////			else (row=idx[j],col=idx[i]);
////			ans+=x[ ((col-1)*N-((col*col-col)>>1)+row-col) -1];
//				// column head index (R) is ((col-1)*N-((col*col-col)>>1)+1
//				// row adjustment is  row-col-1
//				// C indexing adjustment is -1
//				// Note that >> has lower precedence than +/-
//			idx[i]>idx[j] ? (row=idx[i],col=idx[j]) : (row=idx[j],col=idx[i]);
//			ans+=x[ ((col*(N2-col-1))>>1)-N+row-1 ]; // equivalent! but faster?
//		}
//	}
//	ans*=2;
//	return ans;
//}
*/

 double sumSubMatSorted( double const * const x,  int const * const idx,  const int n,  const int N)
/*
// Input:
//   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
//	 idx is a vector of "sorted" indices where sum of x is taken over (idx starts from 1 instead of 0);
//   n is the length of idx; 
//   *N is the total sample size; 
// Output: 
//   returns sum(as.matrix(x)[idx, idx])
// Note:
//   This is about 3~5% faster than sumSubMat when N=30 and B=5000. 
*/
{
	double ans=0.0;
	register unsigned int j, i, N2=(N<<1);
	register int base;

	for(j=0; j<n-1; ++j){ /* // column index */
		/*  base=((idx[j]*(N2-idx[j]-1))>>1)-N-1;  /* // part that does not involve row index */
		base = (((N2-idx[j])*(idx[j]-1))>>1) - idx[j] - 1 ;  /* this should replace the previous line  */
		for(i=j+1; i<n; ++i){  /* //  row index */
			ans+= x [base + idx[i]]; 
		}
	}
/*	ans*=2.0;	*/
	return ans*2.0;
}

#ifdef DEBUG
void testSumSubMatSorted(double const * const x, int const * const idx, const int * const n, const int * const N, double  * const ans)
{
	*ans=sumSubMatSorted(x, idx, *n, *N);
}
#endif

void mrppstats2(double const * const x, int const * const perm, int const * const permcomplement, 
			   int const * const n, int const * const B, int const * const N, int const * const wtmethod,
			   double * const ans)
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
{	register unsigned int b, Nn;
	double tmp1, tmp2, denom1, denom2;

	Nn = *N - *n;
	denom1=1.0/(*n-*wtmethod); denom2=1.0/(Nn-*wtmethod) ;
	for(b=0; b<*B; ++b) 
	{
/*
//		tmp1=sumSubMat(x, perm+(b**n), *n, *N); 
//		tmp2=sumSubMat(x, permcomplement+(b*(*N-*n)), *N-*n, *N); 
*/
		tmp1=sumSubMatSorted(x, perm+(b**n), *n, *N); 
		tmp2=sumSubMatSorted(x, permcomplement+(b*Nn), Nn, *N); 
		

   	    ans[b]=tmp1*denom1 +  tmp2*denom2 ;
	}
}
