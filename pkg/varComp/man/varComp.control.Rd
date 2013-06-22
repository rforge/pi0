\name{varComp.control}
\alias{varComp.control}
\alias{nlminb.control}
\alias{informationTypes}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions controlling the model fitting of varComp
}
\description{
\code{varComp.control} provides detailed control of model fitting for \code{\link{varComp}}. \code{nlminb.control} provides detailed control of the optimization procedure. \code{informationTypes} is a character vector containing the types of information matrices.
}
\usage{
varComp.control(verbose = FALSE, start = NULL, REML = TRUE, 
               information = informationTypes, boundary.eps = 5e-04, 
			   nlminb = nlminb.control(iter.max = 200L, eval.max = 500L), 
			   plot.it = FALSE, keepXYK = TRUE)
nlminb.control(eval.max = 200L, iter.max = 150L, trace = 0L, abs.tol = 0, 
               rel.tol = 1e-10, x.tol = 1.5e-08, xf.tol = 2.2e-14, step.min = 1, 
			   step.max = 1, sing.tol = rel.tol, ...)
informationTypes
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{verbose}{
	A logical scalar, indicating whether excessive messages are printed. 
}
  \item{start}{
	An optional numeric vector of starting values for the ratio of variance components to the error variance. If missing, \code{\link{minque}} will be called with zero as the prior value and with non-negativity constraints. The result of \code{\link{minque}} will be used as the starting value. 
}
  \item{REML}{
	A logical scalar. Currently, only \code{TRUE} is supported. 
}
  \item{information}{
	A character vector (only the first element will be used) specifying the information matrix to be used during optimization. Available choices are given in \code{informationTypes}. See details. 
}
  \item{boundary.eps}{
	A small positive number. If the parameter estimates fall below this number, additional checking will be performed to examine if the parameter estimate hits the lower boundary (0) of the parameter space. 
}
  \item{nlminb}{
	A list of arguments to be passed to \code{\link[stats:nlminb]{stats::nlminb}}. This is most conveniently specified through \code{\link{nlminb.control}}. 
}
  \item{plot.it}{
	A logical scalar, indicating whether the profiled REML function is plotted around the fitted maximum. This is only effective when the number of variance component is one (not counting the error variance). 
}
  \item{keepXYK}{
	A logical scalar, indicating whether the fixed-effect design matrix (\eqn{\mathbf{X}}{X}), response variable (\eqn{\mathbf{Y}}{Y}) and the list of variance-covariance matrices (\eqn{\mathbf{K}}{K}) are included in the result. See \code{\link{varComp}} for notations. 
}  
  \item{eval.max, iter.max, trace, abs.tol, rel.tol, x.tol, xf.tol, step.min, step.max, sing.tol, \dots}{ The same as in \code{\link[stats:nlminb]{stats::nlminb}}. 
}
}
\details{
The objective function is always the profiled restricted log likelihood function, where the error variance is profiled out from the REML function. The parameterization being used is the ratio of each variance component to the error variance. Non-negativity constraints are imposed on such ratios. 

The \code{information} argument specifies the (approximate) Hessian matrix to be used during optimization. Available choices are \code{"EI"} (expected information), \code{"OI"} (observed information), \code{"AEI"} (approximate expected information), \code{"AOI"} (approximate observed information), and \code{"WAI"} (weighted average information). The approximate versions of \code{"EI"} and \code{"OI"} replace the more computationally demanding trace terms by their realized counter-parts that are unbiased for the trace, and the \code{"WAI"} is similar to the average-information of Gilmour (1995) that cancels such terms, except that the latter reference is based on unprofiled REML function rather than profiled REML. The choice of information only affects the speed and the stability of fitting. As long as the objective function is well behaved, the parameter estimates should not be affected by which information matrix is used during fitting. 

Both \code{varComp.control} and \code{nlminb.control} are simply a wrapper of the input arguments. 
}
\value{
A named list of the input arguments. 
}
\references{
Arthur R. Gilmour, Robin Thompson, Brian R. Cullis (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models.  \emph{Biometrisc} 51, 1440--1450. 

	Qu L, Guennel T, Marshall SL. (2013+) Linear Score Tests for Variance Components in Linear Mixed Models and Applications to Genetic Association Studies. \emph{Biometrics}, under revision. 
}
\author{
Long Qu
}
\note{
Setting \code{iter.max=0L} is a trick to get residual contrast and working correlation matrices computed, but not fitting the models completely. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{varComp}}
}
\examples{

library(nlme)
data(Oxide)
## Using fisher scoring algorithm:
varComp(Thickness~Source, Oxide, ~Lot/Wafer, 
  control=varComp.control(information='EI')) 
## Using Newton-Raphson algorithm: 
varComp(Thickness~Source, Oxide, ~Lot/Wafer, 
  control=varComp.control(information='OI')) 
## Computing MINQUE0, residual contrasts, and marginal correlations of residual contrasts
varComp(Thickness~Source, Oxide, ~Lot/Wafer, 
  control=varComp.control(nlminb=nlminb.control(iter.max=0L)))
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize}% __ONLY ONE__ keyword per line
