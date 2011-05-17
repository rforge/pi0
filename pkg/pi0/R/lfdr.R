lfdr=ppee=function(object)UseMethod('lfdr')

lfdr.parncpt=ppee.parncpt=
lfdr.sparncpt=ppee.sparncpt=
lfdr.nparncpt=ppee.nparncpt=
function(object)
{
    pmin(pmax(object$pi0*dt(object$data$tstat, object$data$df)/fitted(object), 0), 1)
}
