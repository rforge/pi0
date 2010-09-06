
//void unisqeuc(double *x, int *nr,  double *dist)
/*********** CALL FROM R ***************
.C("unisqeuc", x = as.double(1:2326), nr = as.integer(2326), dist=dist.holder, DUP = FALSE,  NAOK = TRUE)$dist 
	where x is univariate double vector; 
	     nr is a scalar integer equal to the length x; 
		 dist is a double vector of length nr*nr, which would store the resulting squared euclidean distance matrix
***************************************/
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

inline double sumSubMat(double *x, int *idx, int n, int N)
// Input:
//   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
//	 idx is a vector of sorted indices where sum of x is taken over;
//   n is the length of idx; 
//   *N is the total sample size; 
// Output: 
//   returns sum(as.matrix(x)[idx, idx])
{
	double ans=0.0;
	int i,j,row,col;

	for(i=1; i<n; ++i){
		for(j=0; j<i; ++j){
			if(idx[i]>idx[j]) row=idx[i],col=idx[j];
			else row=idx[j],col=idx[i];

			ans+=x[ ((col-1)*N-((col*col-col)>>1)+row-col) -1];
		}
	}
	ans*=2;
	return ans;
}
void mrppstats(double *x, int *perm, int *permcomplement, int *n, int *B, int *N, int *wtmethod,double *ans)
// Input:
//   x is a stacked vector of the lower triangle of distance mat of (*N x *N); 
//   perm is perm idx mat for 1st group (*n x *B); 
//   permcomplement is perm idx mat for 2nd group ((*N-*n) x *B); 
//   *n is the number of rows in perm; *B is the number of cols in perm; 
//   *N is the total sample size; 
//   *wtmethod is the treatment group weight: 0=sample size-1; 1=sample size
// Output: 
//   *ans is a length *B space of mrpp statistics to be filled in
{	int b;
	double tmp, tmp1;
	for(b=0; b<*B; ++b) {
		tmp=sumSubMat(x, perm+(b**n), *n, *N); 
		tmp1=sumSubMat(x, permcomplement+(b*(*N-*n)), *N-*n, *N); 
		
//		printf("sum(d[perm,perm])=%f\tsum(d[-perm,-perm])=%f\n",tmp,tmp1);
		
   	    ans[b]=tmp/(*n-*wtmethod)+tmp1/(*N-*n-*wtmethod) ;
	}
}
