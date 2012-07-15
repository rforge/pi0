mrpp.test.dist <-
function(y, trt, B=choose(length(trt),table(trt)[1]), perm.mat, wtmethod=0, cperm.mat, ...) ## this uses C code
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
    if(missing(cperm.mat)) cperm.mat=apply(perm.mat, 2, function(kk)(1:N)[-kk])
    stats=.C('mrppstats2',y,perm.mat,cperm.mat,nrow(perm.mat),B,N,as.integer(wtmethod[1]), ans=double(B),
            PACKAGE='MRPP',DUP=FALSE)$ans
    ans=list(statistic=c("MRPP statistic"=stats[1]), all.statistics=stats, 
             p.value=mean(stats[1]-stats>=-min(c(1e-8,.5/B))), parameter=c("number of permutations"=B, 'weight method'=wtmethod[1]),
             data.name=dname, .Random.seed=attr(perm.mat,'.Random.seed'),
             method='2-sample MRPP test')
    class(ans)='htest'
    ans
}


mrpp.test.default <-
function(y, ...) {
    ans=mrpp.test.dist(dist(y),...)
    repl.text=paste("dist(", deparse(substitute(y)), ")",sep='')
    ans$data.name=gsub("dist(y)", repl.text, ans$data.name, fixed=TRUE)
    ans
}

mrpp.test.formula <-
function(y, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    mf <- model.frame(y)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
        stop("grouping factor must have exactly 2 levels")
    ans=mrpp.test.dist(dist(mf[[response]]),g,...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " by "))
    ans
}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

