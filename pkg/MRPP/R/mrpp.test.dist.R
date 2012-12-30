mrpp.test.dist <-
function(y, trt, B=as.integer(min(nparts(table(trt)), 1e4L)), permutedTrt, wtmethod=0, eps=1e-8, ...) ## this uses C code
## y is a dist object; wtmethod: 0=sample size-1; 1=sample size
{
    if(missing(y) || !inherits(y,'dist')) stop('dist object missing or incorrect')
    N= as.integer( ( 1 + sqrt(1.0 + 8.0 * length(y))) * .5   +.5);          ### i.e.,   attr(y,'Size'), howerver, attributes might be lost during subsetting. 
    if(missing(trt)) {  ## recoving trt from the first permutation
      trt=trt.permutedTrt(permutedTrt)
    }
    if(missing(permutedTrt)) {
        permutedTrt=permuteTrt(trt,B, ...)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
    B=nperms.permutedTrt(permutedTrt)
    #if(missing(cperm.mat)) cperm.mat=apply(permutedTrt, 2, function(kk)(1:N)[-kk])
    tabtrt=table(trt)
    ntrt=length(tabtrt)
    stats=.Call(mrppstats,y,permutedTrt, as.numeric(wtmethod[1L]), PACKAGE='MRPP')
#    stats0=.Call(mrppstats, y, 
#        lapply(names(tabtrt), function(z)matrix(which(as.character(trt)==z))), 
#        as.numeric(wtmethod[1L]),         PACKAGE='MRPP')
    ans=list(statistic=c("MRPP statistic"=stats[1]), all.statistics=stats, 
             p.value=mean(stats[1]-stats>=-eps), parameter=c("number of permutations"=B, 'weight method'=wtmethod[1L]),
             data.name=dname, #  .Random.seed=attr(permutedTrt,'.Random.seed'),  ## this is commented out since the random number seed in gmp:::urand.bigz will be used. 
             method=sprintf('%d-sample MRPP test',ntrt)
             )
    class(ans)='htest'
    ans
}


mrpp.test.default <-
function(y, ...) {
    ans=mrpp.test.dist(dist(y),...)
    repl.text=paste("Response data ", deparse(substitute(y)), sep='')
    ans$data.name=gsub("\"dist\" object dist(y)", repl.text, ans$data.name, fixed=TRUE)
    ans
}

mrpp.test.formula <-
function(y, data, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    if(missing(data)) mf = model.frame(y) else mf <- model.frame(y, data=data)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    res=mf[[response]]
    ans=mrpp.test(dist(res),trt=g, ...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " ~ "))
    if(!missing(data)) ans$data.name = paste(ans$data.name, "in data.frame", substitute(data))
    ans
}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

