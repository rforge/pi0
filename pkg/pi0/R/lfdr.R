lfdr=ppee=function(object, ...)UseMethod('lfdr')

lfdr.parncpt=ppee.parncpt=
lfdr.sparncpt=ppee.sparncpt=
lfdr.nparncpt=ppee.nparncpt=
function(object, ...)
{
    pmin(pmax(object$pi0*dt(object$data$tstat, object$data$df)/fitted(object), 0), 1)
}

lfdr.nparncpp=ppee.nparncpp=
function(object, ...)
{
    object$LFDR
}

lfdr.CBUM=ppee.CBUM=
function(object, ...)
{
    attr(object, 'lfdr')
}

lfdr.znormix=ppee.znormix=
function(object, ...)
{
    attr(object, 'lfdr')
}

lfdr.convest=ppee.convest=
function(object, ...)
{
    attr(object, 'lfdr')
}
