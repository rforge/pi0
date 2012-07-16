mrpp.test.dist <-
function(y, trt, B=choose(length(trt),table(trt)[1]), perm.mat, wtmethod=0, eps=1e-8, ...) ## this uses C code
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
    #if(missing(cperm.mat)) cperm.mat=apply(perm.mat, 2, function(kk)(1:N)[-kk])
    tabtrt=table(trt)
    ntrt=length(tabtrt)
    stats=.Call('mrppstats',y,perm.mat, as.numeric(wtmethod[1L]), PACKAGE='MRPP')
    stats0=.Call('mrppstats', y, 
        lapply(names(tabtrt), function(z)matrix(which(as.character(trt)==z))), 
        as.numeric(wtmethod[1L]),         PACKAGE='MRPP')
    ans=list(statistic=c("MRPP statistic"=stats0), all.statistics=stats, 
             p.value=mean(stats0-stats>=-eps), parameter=c("number of permutations"=B, 'weight method'=wtmethod[1L]),
             data.name=dname, .Random.seed=attr(perm.mat,'.Random.seed'),
             method=sprintf('%d-sample MRPP test',ntrt)
             )
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
function(y, data, B, perm.mat, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    if(missing(data)) mf = model.frame(y) else mf <- model.frame(y, data=data)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    resNames=mf[[response]]
    res=do.call('cbind', lapply(resNames, 'get', pos=if(missing(data)) -1 else data))
    ans=mrpp.test.dist(dist(res),g,...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " ~ "))
    if(!missing(data)) ans$data.name = paste(ans$data.name, "in data.frame", substitute(data))
    ans
}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

