\name{   ncpest       }
\alias{    ncpest     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
semipar NCP estimate
}
\description{
semipar NCP estimate

}
\usage{
ncpest(p,
    breaks=min(2000,round(length(p)/5)),
    test=c("t","z"),
    df,
    alternative=c("two.sided", "less", "greater"),
    compromise.n=1,
    lambdas=#if(penalty_type==1)10^seq(-2,6,length=6) else 
            10^seq(-4,6,length=11),
    deltamax='auto',
    nknots,
    ndelta=500,
    solver=c("LowRankQP","lsei","solve.QP","ipop"),
    weights=1,
    keep.cdf=TRUE,
    LowRankQP.method=c('LU','CHOL'),
    lsei.method=c('chol','svd','eigen'),
    debugging=FALSE,
    ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{              }
  \item{breaks}{              }
  \item{test}{              }
  \item{df}{              }
  \item{alternative}{              }
  \item{compromise.n}{              }
  \item{lambdas}{              }
  \item{deltamax}{              }
  \item{nknots}{              }
  \item{ndelta}{              }
  \item{solver}{              }
  \item{weights}{              }
  \item{keep.cdf}{              }
  \item{LowRankQP.method}{              }
  \item{lsei.method}{              }
  \item{debugging}{              }
  \item{\dots}{              }
}
\details{





}
\value{
%  \item{}{              }
%  \item{}{              }
%  \item{}{              }
  
}
\references{


}
\author{ Long Qu [\email{longor@iastate.edu}] }
\note{


}
\seealso{ 


}
\examples{
\dontrun{


}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar  }
%\keyword{  }
%\keyword{  }
