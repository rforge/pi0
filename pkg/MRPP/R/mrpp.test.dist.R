mrpp.test.dist <-
function(y, trt, B=choose(length(trt),table(trt)[1]), perm.mat) ## this is C code
## y is a dist object
{
    if(missing(y) || !inherits(y,'dist')) stop('dist object missing or incorrect')
    N=attr(y,'Size')
    if(missing(perm.mat)) {
        perm.mat=get.perm.mat(trt,B)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permutation matrix', deparse(substitute(perm.mat)))
    B=ncol(perm.mat)
    cperm.mat=apply(perm.mat, 2, function(kk)(1:N)[-kk])
    stats=.C('mrppstats',y,perm.mat,cperm.mat,nrow(perm.mat),B,N,ans=double(B))$ans
    ans=list(statistic=c("MRPP statistic"=stats[1]), all.statistics=stats, 
             p.value=mean(stats[1]>=stats), parameter=c("number of permutations"=B),
             data.name=dname, .Random.seed=attr(perm.mat,'.Random.seed'),
             method='2-sample MRPP test')
    class(ans)='htest'
    ans
}

