\name{   nparncpp       }
\alias{    nparncpp     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Estimation of the density of absolute noncentrality parameters
}
\description{
Estimation of the density of absolute noncentrality parameters, using linear B-spline model. 
}
\usage{
nparncpp(p,
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
    solver=c("lsei","LowRankQP","solve.QP","ipop"),
    weights=1,
    keep.cdf=TRUE,
    LowRankQP.method=c('LU','CHOL'),
    lsei.method=c('chol','svd','eigen'),
    debugging=FALSE,
    ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{p-value vector}
  \item{breaks}{break points to bin the p-values}
  \item{test}{either \code{t}-test or \code{z}-test}
  \item{df}{degrees of freedom for the test}
  \item{alternative}{Same as in \code{\link{t.test}}}
  \item{compromise.n}{Number of components in the compromised estimate}
  \item{lambdas}{Candidate tuning parameters}
  \item{deltamax}{Assumed maximum noncentrality parameters}
  \item{nknots}{Number of knots}
  \item{ndelta}{Number of points to evaluate the noncentrality parameters}
  \item{solver}{Quadratic programming solver function}
  \item{weights}{Bin weights}
  \item{keep.cdf}{Logical: whether computed conditional CDF is saved in global enviroment. See \code{\link{cond.cdf}}. }
  \item{LowRankQP.method}{Method for \code{\link[LowRankQP:LowRankQP]{LowRankQP}}}
  \item{lsei.method}{Method for \code{\link[limSolve:lsei]{lsei}}}
  \item{debugging}{Logical: print excessive messages}
  \item{\dots}{Additional argumenets to \code{solver}}
}
%\details{
%
%
%
%
%
%}
\value{An  object of class \code{c('nparncpp','ncpest')}.
%  \item{}{              }
%  \item{}{              }
%  \item{}{              }
  
}
\references{
Ruppert, Nettleton, Hwang. (2007) Exploring the Information in \eqn{p}-values for the Analysis and Planning of Multiple-test Experiments. Biometrics. 63. 483-495.

}
\author{ Long Qu  translated and enhanced the original \code{MATLAB} code from Dr. David Ruppert.}
\note{
The code right now is not completely compatible with the \code{ncpest} class and is subject to change in future versions.

}
\seealso{ \code{\link{nparncpt}}, \code{\link{sparncpt}}, \code{\link{parncpt}}, \code{\link{dncp}} 

}
%\examples{
%\dontrun{
%
%
%}
%
%
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth  }
%\keyword{  }
%\keyword{  }
