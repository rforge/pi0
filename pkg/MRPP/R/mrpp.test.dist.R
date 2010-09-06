mrpp.test.dist <-
function(y, trt, B=choose(length(trt),table(trt)[1]), perm.mat, wtmethod=0) ## this is C code
## y is a dist object; wtmethod: 0=sample size-1; 1=sample size
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
    stats=.C('mrppstats',y,perm.mat,cperm.mat,nrow(perm.mat),B,N,as.integer(wtmethod[1]), ans=double(B))$ans
    ans=list(statistic=c("MRPP statistic"=stats[1]), all.statistics=stats, 
             p.value=mean(stats[1]>=stats), parameter=c("number of permutations"=B, 'weight method'=wtmethod[1]),
             data.name=dname, .Random.seed=attr(perm.mat,'.Random.seed'),
             method='2-sample MRPP test')
    class(ans)='htest'
    ans
}

mrpp.test.matrix <-
function(y, ...) {
    if(deparse(substitute(y))=='y'){
        mrpp.test.dist(as.dist(y),...)
    }else{
        assign(deparse(substitute(y)),y)
        eval(substitute(mrpp.test.dist(as.dist(yyy),...), list(yyy=substitute(y))))
    }
}

mrpp.test.formula <-
function(y,B, perm.mat) ## not implemented
{

}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

