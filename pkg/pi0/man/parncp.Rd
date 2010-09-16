\name{parncpt}
\Rdversion{1.1}
\alias{parncpt}
\alias{parncpt.bfgs.0mean}
\alias{parncpt.bfgs.non0mean}
\alias{parncpt.momeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric estimation of noncentrality parameter distribution
}
\description{Assuming normality of noncentrality parameters, the MLE of its standard deviation (and possibly mean also) is estimated from observed t-statistics
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parncpt(tstat, df, zeromean = TRUE, ...)
parncpt.bfgs.0mean(tstat, df, starts, grids, approximation = "int2", ...)
parncpt.bfgs.non0mean(tstat, df, starts, grids, approximation = "int2", ...)
parncpt.momeff(tstat,n1,n2=n1,zeromean,gamma2,lower.df=6.1,upper.df=100,approx=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{numeric vector of t-statistics
%%     ~~Describe \code{tstat} here~~
}
  \item{df}{numeric vector of degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{zeromean}{logical; if \code{TRUE}, then mean of noncentrality parameters is assumed to be zero and is \emph{not} estimated. 
%%     ~~Describe \code{zeromean} here~~
}
  \item{\dots}{Other arguments to \code{\link{optim}}
%%     ~~Describe \code{\dots} here~~
}
  \item{starts}{An optional vector of starting values. If missing, a grid search will be performed to get a good starting value. 
%%     ~~Describe \code{starts} here~~
}
  \item{grids}{A list of three components (\code{lower}, \code{upper}, \code{ngrid}) defining the grids to be searched in find a good starting value. 
  Each component is a numeric vector of the same length as the number of parameters. \code{lower} and \code{upper} give the bounds, and
  \code{ngrid} specifies the number of points for each dimension. 
%%     ~~Describe \code{grids} here~~
}
  \item{approximation}{Methods of approximating the noncentral t-density. \code{int2} is exact for integer df, but interpolate to fractional df. 
  'laplace' is the laplacian approximation; 'saddlepoint' is the saddlepoint approximation; 'none' computes the (sort of) exact density using the default \code{\link{dt}} function. 
%%     ~~Describe \code{approximation} here~~
}
  \item{n1}{Treatment 1 sample size}
  \item{n2}{Treatment 2 sample size}
  \item{gamma2}{Gamma square parameter, i.e., variance of effect sizes.}
  \item{lower.df}{  lower bound of degrees of freedom, in case of n1 is missing            }
  \item{upper.df}{  upper bound of degrees of freedom, in case of n1 is missing            }
  \item{approx}{  logical, indicating if no exact solutions are available, whether approx. solutions are returned.  }

}
\details{\code{parncpt} calls either \code{parncpt.bfgs.0mean} or \code{parncpt.bfgs.non0mean}, depending whether \code{zeromean} is \code{TRUE} or \code{FALSE}. 
Both \code{parncpt.bfgs.0mean} and \code{parncpt.bfgs.non0mean} use the 'L-BFGS-B' algorithm by calling \code{\link{optim}}. All gradiants are analytical, but the Hessian is only numerical approximation. 
The first parmater is always \code{pi0}, i.e., the proportion of true null hypotheses; the last parameter is always the standard deviation of noncentrality parameters; 
for \code{parncpt.bfgs.non0mean} the middle parameter is the mean of noncentrality parameters, whereas for \code{parncpt.bfgs.0mean} the mean is set to 0 a priori. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{a list with \code{class} attribute being \code{c('parncpt', 'ncpest')}. 
        \item{pi0}{proportion of true nulls}
        \item{mu.ncp}{mean of ncp}
        \item{sd.ncp}{SD of ncp}
        \item{data}{a list of \code{tstat} and \code{df}}
        \item{logLik}{an object of class \code{logLik}. Call \code{\link{logLik.ncpest}} to extract. Similarly, \code{\link{AIC}} is callable.}
        \item{enp}{the (effective) number of parameters in the model}
        \item{par}{estimated parameters. Call \code{\link{coef.ncpest}} to extract.}
        \item{obj}{the negative loglikelihood function that is minimized}
        \item{gradiant}{analytic gradiant at the estimate}
        \item{hessian}{numeric hessian at the estimate}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Long Qu, Dan Nettleton, Jack CM Dekkers. Improved estimation of noncentrality parameter distribution from a large number of t-statistics, with applications to false discovery rates in microarray data. (manuscripts under preparation). 
%% ~put references to the literature/web site here ~
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
\note{\code{df} could be \code{Inf} for z-tests. When this is the case, \code{approximation} is ignored. 

\code{parncpt.momeff} is the old code using method of moments estimates. It is outdated, depreciated, and not completely compatible with current \code{ncpest} class. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sparncpt}}, \code{\link{nparncpt}}, 
\code{\link{fitted.parncpt}}, \code{\link{plot.parncpt}}, \code{\link{summary.parncpt}},
\code{\link{coef.ncpest}}, \code{\link{logLik.ncpest}}, \code{\link{vcov.ncpest}},
\code{\link{AIC}}, \code{\link{dncp}}
}
\examples{
\dontrun{
data(simulatedTstat)
(npfit=nparncpt(tstat=simulatedTstat, df=8)); 
(pfit=parncpt(tstat=simulatedTstat, df=8, zeromean=FALSE)); plot(pfit)
(pfit0=parncpt(tstat=simulatedTstat, df=8, zeromean=TRUE)); plot(pfit0)
(spfit=sparncpt(npfit,pfit)); plot(spfit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
