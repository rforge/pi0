#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void tstatistic(double *x,int *n1,int *n2,int *ntests,int *byrow,int *pool,double *tstat, double *df)
{
     double mean1,mean2,S1,S2,pooledvar;
     int i,l,ncol,n=(*n1)+(*n2);
/*
	if(((*n1))<2 || *n2<2){
              MessageBox(0,"n is too small","ERROR",1);
              exit(0);
     }
*/

     if (*byrow==1){
       for(l=0;l<*ntests;l++){
         mean1=mean2=S1=S2=pooledvar=0.0;

         for(i=0;i<(*n1);i++)mean1+=x[l+i**ntests];
         mean1/=(double)(*n1);                           
         for(i=(*n1);i<n;i++)mean2+=x[l+i**ntests];
         mean2/=(double)(*n2);                           
         
         S1=S2=0.0;
         for(i=0;i<(*n1);i++){
             S1+=(x[l+i**ntests]-mean1)*(x[l+i**ntests]-mean1);
         }
         for(i=(*n1);i<n;i++){
             S2+=(x[l+i**ntests]-mean2)*(x[l+i**ntests]-mean2);
         }

         if(*pool){
                   pooledvar=(S1+S2)/(n-2);
                   tstat[l]=(mean1-mean2)/
                     sqrt( pooledvar * (1/(double)((*n1))+1/(double)(*n2) )  );
         }else{
               tstat[l]=(mean1-mean2)/
                     sqrt(S1/(double)((*n1)*((*n1)-1))+S2/(double)(*n2*((*n2)-1)));
                S1/=(double)((*n1)-1);
                S2/=(double)((*n2)-1);
               df[l]=pow(S1/(double)(*n1)+S2/(double)(*n2),2.0) /
                    (S1*S1/(double)((*n1)*(*n1)*((*n1)-1)) + 
                     S2*S2/(double)((*n2)*(*n2)*((*n2)-1)));
         } /*of pooling variance*/
       } /*of one line*/
     } /*of byrow*/
     else{
       ncol=*ntests;  
       for(l=0;l<ncol;l++){
         mean1=mean2=S1=S2=pooledvar=0.0;

         for(i=0;i<(*n1);i++)mean1+=x[l*n+i];
         mean1/=(double)(*n1);
         for(i=(*n1);i<n;i++)mean2+=x[l*n+i];
         mean2/=(double)(*n2);
         
         for(i=0;i<(*n1);i++)
             S1+=(x[l*n+i]-mean1)*(x[l*n+i]-mean1);
         for(i=(*n1);i<n;i++)
             S2+=(x[l*n+i]-mean2)*(x[l*n+i]-mean2);

         if(*pool){
                   pooledvar=(S1+S2)/(n-2);
                   tstat[l]=(mean1-mean2)/
                     sqrt( pooledvar * (1/(double)((*n1))+1/(double)(*n2) )  );
         }else{
               tstat[l]=(mean1-mean2)/
                     sqrt(S1/(double)((*n1)*((*n1)-1))+S2/(double)(*n2*(*n2-1)));
                S1/=(double)((*n1)-1);
                S2/=(double)((*n2)-1);
               df[l]=pow(S1/(double)(*n1)+S2/(double)(*n2),2.0) /
                    (S1*S1/(double)((*n1)*(*n1)*((*n1)-1))
                    +S2*S2/(double)((*n2)*(*n2)*((*n2)-1)));
         } /*of pooling variance*/
       } /*of one line*/
     } /*of byrow else*/
 
}


