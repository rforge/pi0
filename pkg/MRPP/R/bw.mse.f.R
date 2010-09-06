bw.mse.F <-
function(z, mc.n=500)
{
    K=mc.n
    y=qnorm(1:K/(K+1))
    Phi.y=pnorm(y)
    t=z[1]
    h0=bw.nrd(z)
    B=length(z)
    t.z.h0=(t-z)/h0
    y.h0=y/h0
    B.F.0=sum(pnorm(t.z.h0))
    inner=Phi.inner=sum.Phi.inner=current.h=Inf    ## to be modified by obj or der.obj
    obj=function(h){
        if(current.h!=h){
            current.h<<-h
            inner<<-outer(t.z.h0, h*y.h0, '-')
            Phi.inner<<-pnorm(inner)
            sum.Phi.inner<<-sum(Phi.inner)
        }
      2/B/B/K*sum(sweep(Phi.inner, 2, Phi.y, '*'))+
      (B-1)/B/B/B/K/K*sum.Phi.inner*sum.Phi.inner-
      2/B/B/K*B.F.0*sum.Phi.inner
    }
    der.obj=function(h){
        if(current.h!=h){
            current.h<<-h
            inner<<-outer(t.z.h0, h*y.h0, '-')
            Phi.inner<<-pnorm(inner)        # not used in der.obj, but keep them up-to-date
            sum.Phi.inner<<-sum(Phi.inner)  # not used in der.obj, but keep them up-to-date
        }
        phi.inner=dnorm(inner)
        y.phi.inner=sweep(phi.inner,2,y,'*')
            sum.y.phi.inner=sum(y.phi.inner)
        Phi.y.phi.inner=sweep(y.phi.inner,2,Phi.y,'*')
     -2/B/B/h0/K*sum(Phi.y.phi.inner)-
      2*(B-1)/B/B/B/h0/K/K*sum.Phi.inner*sum.y.phi.inner+
      2/B/B/h0/K*B.F.0*sum.y.phi.inner
    }
    ofit=optim(bw.mse.f.asym(z,iter=FALSE), obj, der.obj, method='L-BFGS-B', lower=0, 
                control=list(fnscale=1e-6#,factr=1e3,parscale=1e-5
                ))
    ans=ofit$par
    attr(ans,'optim')=ofit
    ans
}

