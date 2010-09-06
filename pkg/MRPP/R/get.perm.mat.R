get.perm.mat <-
function(trt, B=100) ## trt needs have exactly 2 levels; other situations are not implemented yet
{
    if(!exists('.Random.seed', envir=globalenv())) runif(1)
    save.seed=get('.Random.seed', envir=globalenv())
    n=table(trt)
    perms=replicate(B, sort(sample(sum(n),min(n))))
    perms[,1]=which(as.factor(trt)==names(n[which.min(n)]))
    attr(perms,'.Random.seed')=save.seed
    perms
}

