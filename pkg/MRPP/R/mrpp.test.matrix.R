mrpp.test.matrix <-
function(y, ...) {
    if(deparse(substitute(y))=='y'){
        mrpp.test.dist(as.dist(y),...)
    }else{
        assign(deparse(substitute(y)),y)
        eval(substitute(mrpp.test.dist(as.dist(yyy),...), list(yyy=substitute(y))))
    }
}

