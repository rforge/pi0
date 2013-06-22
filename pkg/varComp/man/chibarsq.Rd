\name{Chi-bar-square}
\alias{mchibarsq}
\alias{pchibarsq}
\alias{wchibarsq}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Chi-bar-square distribution with nonnegativity cone constraints
}
\description{
\code{pchibarsq} is the distribution function of chi-bar-square distribution with nonnegativity cone constraint. 

\code{wchibarsq} computes the mixing proportions for the chi-bar-square distribution.

\code{mchibarsq} computes the moments of the chi-bar-square distribution. 
}
\usage{
pchibarsq(q, V, lower.tail = TRUE, log.p = FALSE)
mchibarsq(V, order = 1:2)
wchibarsq(V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{q}{ A vector of quantiles, as in \code{\link[stats:pchisq]{stats::pchisq}}.
}
  \item{V}{ A positive-definite matrix, defining the distance measure used when projecting onto the cone. 
}
	\item{lower.tail}{logical, the same as in \code{\link[stats:pchisq]{stats::pchisq}}. 
}
	\item{log.p}{logical, the same as in \code{\link[stats:pchisq]{stats::pchisq}}. 
}
  \item{order}{
	A positive integer vector of the order of moments to be computed.
}
}
% \details{
% }
\value{
\code{pchibarsq} gives the distribution function, \code{wchibarsq} gives the mixing proportion,  and \code{mchibarsq} gives the moments, 
}
\references{
A. Shapiro (1988) Towards a Unified Theory of Inequality Constrained Testing in Multivariate Analysis. \emph{Int. Stat. Rev.} 56, 49--62. 

Akio Kudo (1963) A multivariate analogue of the one-sided test.
 \emph{Biometrika} 50, 403--418.
}
\author{
	Long Qu
}

\seealso{
	\code{\link[stats:pchisq]{stats::pchisq}}
}
\examples{
set.seed(203490L)
  V=crossprod(matrix(rnorm(25),5))
  VI=solve(V)
  L=t(chol(V))
  chibarsq=replicate(1e3L, -2*quadprog::solve.QP(VI, VI\%*\%(L\%*\%rnorm(5)),
     diag(1,5), rep(0,5))[['value']] )
  chibarsq=sort(chibarsq)
  p=pchibarsq(chibarsq, V)
  
\dontrun{
  plot(ecdf(chibarsq))
  lines(chibarsq, p, col=4, lwd=3, lty=3)
}
  mean(chibarsq); mean(chibarsq^2)
  mchibarsq(V)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
