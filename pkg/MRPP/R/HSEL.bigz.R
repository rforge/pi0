HSEL.bigz=function(N, M, seed=0)
#   Sample M out of N without replacement where N can be a huge integer. 
#    Ref: Jarmo Ernvall and Olli Nevalainen. 1982. An Algorithm for Unbiased Random Sampling. THE COMPUTER JOURNAL, VOL. 25, NO. 1, 45--47
{

    M=as.integer(M)
    stopifnot(M<=N)
    if(M==0L) return(as.bigz(integer(0L)))
    IS=IA=IB=as.bigz(integer(M))
    IH=LINK=integer(M)
    LAST=N
    NEW=0L
    urand.bigz(1, seed=seed)
    for(L in 1:M){
        KSI=sample1(LAST)+1
        I=as.integer(mod.bigz(LAST, M) + 1)
        J=as.integer(IH[I])
        repeat{ # added line 4
            if (J == 0L) break # line 4
            if (IA[J] == LAST) break
            J=LINK[J]
        }
        MLAST=LAST # line 3
        if (J != 0L) MLAST=IB[J]
        I=as.integer(mod.bigz(KSI, M) + 1)
        J=as.integer(IH[I])
        flag7=FALSE; flag8=FALSE
        repeat{
            if (J == 0) {flag7=TRUE; break }# line 9; # to line 7
            if (IA[J] == KSI) {flag8=TRUE; break}
            J=LINK[J]
        }
        if(flag8){
            IS[L]=IB[J]
            IB[J]=MLAST
            flag7=FALSE
        }
        if(flag7){
            IS[L]=KSI   # line 7
            NEW=NEW+1
            IA[NEW]=KSI
            IB[NEW]=MLAST
            LINK[NEW]=IH[I]
            IH[I]=NEW
        }
        LAST = LAST - 1 # line 2
    }
    as.vector(IS)
}

sample1=function(N, seed=0)
# sample one number between 0 and N-1, where N can be huge
{
    expo=frexpZ(N)$exp
    ans = urand.bigz(1, expo, seed)
    while(ans >= N)        ans = urand.bigz(1, expo)
    ans
}
