#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP combnBtRc(SEXP n, SEXP m, SEXP a, SEXP selected, SEXP R, SEXP nexttoselect, SEXP evalfun, 
			SEXP env, SEXP all)
{
	int next2Sel,k,N=INTEGER(n)[0],M=INTEGER(m)[0],*A,*Sel,RR=INTEGER(R)[0],
		goodJump,steps,*next2SelPtr, All=INTEGER(all)[0];
	long nExistAns;
	double jump=0.0,lastJump=0.0;
	A=INTEGER(a);
	Sel=INTEGER(selected);
	next2SelPtr=INTEGER(nexttoselect);

	k=0;
	nExistAns=0;
	next2Sel=0;
	while(1){
/*		printf("(k+1)=%d\tSel[next2Sel]=%d\tnExistAns=%d\n",k+1,Sel[next2Sel],nExistAns);*/
		if(N-M+(k+1)-A[k]<1){
			k--;
			continue;
		}
		if(All){
			A[k]++;
		}else{
			for(jump=0.0,goodJump=0,steps=1;steps<=N-M+(k+1)-A[k];steps++){
				lastJump=jump;
				jump+=choose((double)(N-(A[k]+steps)),(double)(M-(k+1)));
				if(jump>=(Sel[next2Sel]-nExistAns)) {goodJump=1;break;}
			}
			if(!goodJump){		
				nExistAns+=jump;
				k--;			
				continue;		
			}

			A[k]+=steps;
            nExistAns+=lastJump;
		}

        if(k+1==M){
/*			for(i=0;i<M;i++)printf("%d\t",A[i]);printf("\n");*/
            nExistAns++;
			next2SelPtr[0]=next2Sel+1;
            eval(evalfun,env);
            if(++next2Sel>=RR)break;
        }else{
			++k;
			A[k]=A[k-1]; /* A[++k]=A[k-1]; just to avoid warnings. */
        }
    }

	return(R_NilValue);
}
