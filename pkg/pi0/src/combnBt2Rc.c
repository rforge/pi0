#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


SEXP combnBt2Rc(SEXP n, SEXP m, SEXP n2, SEXP m2, SEXP a, SEXP a2, 
				SEXP selected1, SEXP R1, SEXP selected2All, SEXP R2perGrp1,
				SEXP L,SEXP evalfun, SEXP env, SEXP all)

{
	int N1=INTEGER(n)[0],N2=INTEGER(n2)[0],M1=INTEGER(m)[0],M2=INTEGER(m2)[0],*A1,*A2,
		*Sel1,RR1=INTEGER(R1)[0], *Sel2All, *Sel2, *R2p1, *LL,
		All=INTEGER(all)[0],
		goodJump,steps,i, next2Sel1, next2Sel2,k1,k2;
	long nExistAns1,nExistAns2;
	double jump=0.0,lastJump=0.0;
	A1=INTEGER(a);
	A2=INTEGER(a2);
	Sel2=Sel1=INTEGER(selected1); /* Sel2 is redefined later; just to avoid warnings */
	Sel2All=INTEGER(selected2All);
	R2p1=INTEGER(R2perGrp1);
	LL=INTEGER(L);

	k1=0;
	nExistAns1=0;
	next2Sel1=0;
	while(1){
/*		Rprintf("(k+1)=%d\tnExistAns=%d\n",k1+1,nExistAns1);  */
		if(N1-M1+(k1+1)-A1[k1]<1){
			k1--;
			continue;
		}
		if(All){
			A1[k1]++;
		}else{
			for(jump=0.0,goodJump=0,steps=1;steps<=N1-M1+(k1+1)-A1[k1];steps++){
				lastJump=jump;
				jump+=choose((double)(N1-(A1[k1]+steps)),(double)(M1-(k1+1)));
				if(jump>=(Sel1[next2Sel1]-nExistAns1)) {goodJump=1;break;}
			}
			if(!goodJump){		
				nExistAns1+=jump;
				k1--;			
				continue;		
			}
            nExistAns1+=lastJump;
			A1[k1]+=steps;
		}

        if(k1+1==M1) {

            { /*** group 2 blk ***/

				k2=0;
				nExistAns2=0;
				next2Sel2=0;
				for(i=0;i<M2;i++)A2[i]=0;
				if(!All){
					Sel2=Calloc(*R2p1,int);
					for(i=0;i<*R2p1;i++)Sel2[i]=*(Sel2All++);
				}
				while(1){
					if(N2-M2+(k2+1)-A2[k2]<1){
						k2--;
						continue;
					}
					if(All){
						A2[k2]++;
					}else{
						for(jump=0.0,goodJump=0,steps=1;steps<=N2-M2+(k2+1)-A2[k2];steps++){
							lastJump=jump;
							jump+=choose((double)(N2-(A2[k2]+steps)),(double)(M2-(k2+1)));
							if(jump>=(Sel2[next2Sel2]-nExistAns2)) {goodJump=1;break;}
						}
						if(!goodJump){		
							nExistAns2+=jump;
							k2--;			
							continue;		
						}

						A2[k2]+=steps;
						nExistAns2+=lastJump;
					}
				    if(k2+1==M2){
						nExistAns2++;
						LL[0]++;
/*						Rprintf("%d:\n",LL[0]);
//				for(i=0;i<M1;i++)Rprintf("%d ",A1[i]);Rprintf("\t");
//				for(i=0;i<M2;i++)Rprintf("%d ",A2[i]);Rprintf("\n"); */
						eval(evalfun,env);
						if(++next2Sel2>=*R2p1){
							R2p1++;
							if(!All)Free(Sel2);
							break;
						}
					}else{
						++k2;
						A2[k2]=A2[k2-1]; /* A2[(++k2)]=A2[k2-1];==>this gives warnings. */
					}
				}
			} /*** end of group 2 blk ***/
			/* eval(evalfun,env);*/
            nExistAns1++;
            if(++next2Sel1>=RR1)break;
        }else{
			++k1;
			A1[k1]=A1[k1-1];
        }
    }

	return(R_NilValue);
}


