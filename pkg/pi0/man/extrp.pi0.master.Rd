\name{extrp.pi0}
\alias{extrp.pi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extrapolate the Estimates of P-value Density at 1 from Subsamples to Estimate the Proportion of True Nulls }
\description{
This is the second step of the Subsampling-Extrapolation (SubEx) procedure for estimating the proportion 
of \bold{TRUE} null hypotheses, i.e., \code{pi0}, when a large number of two-sample t-tests are simultaneously 
performed. It regresses the p-value density estimates at 1 from subsamples over various subsample sizes and 
extrapolates the curve/plane to infinite sample sizes in each group. This estimated limit is used to estimate
pi0.
}
\usage{
extrp.pi0(dat,slope.constraint=TRUE,
         gamma2.range=2^c(-4,3),
         rate.margin=c(0,0),
         plotit=TRUE
         ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ an object of class \code{\link{subt}}; typically resulting from calling the function \code{\link{subt}}. }
  \item{slope.constraint}{ logical.}
  \item{gamma2.range}{ a numeric vector of length 2. }
  \item{rate.margin}{ a numeric vector of length 2. }
  \item{plotit}{logical}
}
\details{
  Two regression functions may be used, as specified by \code{nparm}. One is assuming the nonzero standardized effect sizes have a marginal distribution
  of zero mean normal distribution with variance \eqn{\gamma^2}{gamma^2}. This regression function has two parameters,\cr
  \deqn{f_1 =(1-\pi_0)\sqrt{\frac{n_1+n_2}{n_1+n_2+n_1 n_2 \gamma^2}}+\pi_0}{f1=(1-pi0)sqrt[(n1+n2)/(n1+n2+n1*n2*gamma^2)]+pi0,}
  where \eqn{f_1}{f1} is the density estimates at 1 for subsamples, \eqn{n_1}{n1} and \eqn{n_2}{n2} are the corresponding 
  subsampling sizes, \eqn{0 \le \pi_0 \le 1}{0<=pi0<=1}, and \eqn{\gamma^2 > 0}{gamma^2>0}. \cr
  The other regression function is more flexible by replacing the square root and the \eqn{1-\pi_0}{1-pi0} term with another
  two parameters:\cr
  \deqn{f_1 =a \left(\frac{n_1+n_2}{n_1+n_2+n_1 n_2 \gamma^2}\right)^c+\pi_0}{f1=a*[(n1+n2)/(n1+n2+n1*n2*gamma^2)]^c+pi0,} 
  subject to additional constraints of \eqn{a>0}{a>0} and \eqn{c>0}{c>0}.\cr

  The \code{FUN} is the name of an optimization function to be used in non-linear least squares under 
  constrains. Currently, two built-in choices are \code{\link[stats]{constrOptim}} and \code{\link[rgenoud]{genoud}}. 
  If one likes to use either of these two, simply set \code{FUN} to their names. In our experience, 
  \code{\link[stats]{constrOptim}} works well and fast enough. If the resulting estimate is spurious, e.g., 
  too small of too large, one can try the genetic optimization \code{\link[rgenoud]{genoud}} in the \pkg{rgenoud} package. 
  If one does not want to use these two choices, the name of a user defined optimization function can 
  be passed to \code{FUN}. In this case, the function must accept at least two arguments, the \code{dat}
  of class \code{subt} and the numeric vector \code{starts}. Other arguments can be supplied with \code{\dots}.
  Also, the function needs to make sure the estimated parameters are within the parameter space. The value
  of the function needs to be a list, with at least one component named \code{par}, which is a numeric 
  vector of estimated parameters. \cr

  \code{starts} is a numeric vector, the first element of which is the starting value for \eqn{\pi_0}{pi0} 
  and the second element of which is the starting value for \eqn{\gamma^2}{gamma^2}. If \code{nparm=4}, 
  then the third element is the starting value for \eqn{a}{a} and the fourth element is the starting value
  for \eqn{c}{c}. Only the first \code{nparm} elements of \code{starts} are used, with the rest ignored. \cr

  It is highly recommended to have \pkg{rgl} package available to display the estimated regression surface 
  and possibly rotate it with the mouse. 

}
\value{
  an object of class \code{extrpi0}, which is a numeric vector of length 1, named \code{"pi0"}, giving the 
  estimated \eqn{\pi_0}{pi0}, with the following \code{attributes}:\cr

  \item{attr(,'fitted.obj') }{a list, which is the object returned by \code{FUN}. }
  \item{attr(,'nparm') }{ the same as the first element of \code{nparm}. }
  \item{attr(,'extrpFUN') }{the same as \code{FUN}. }
  \item{attr(,'start.val') }{the first \code{nparm} elements of \code{starts}. }
  \item{attr(,'subt.data') }{the same as \code{dat}. }
  
}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }
\note{
  \code{FUN='genoud'} may take a longer time. 
}
\seealso{ 
	\code{\link{subt}}, \code{\link{subex}}, \code{\link[stats]{constrOptim}}, \code{\link[rgenoud]{genoud}}, 
	\code{\link{print.extrpi0}}, \code{\link{plot.extrpi0}}. 
}
\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)	
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf)
## this is how the 'simulatedExtrpi0' data set in this package generated
simulatedExtrpi0=extrp.pi0(simulatedSubt)
}
data(simulatedExtrpi0)
summary(simulatedExtrpi0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ htest }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{optimize}
