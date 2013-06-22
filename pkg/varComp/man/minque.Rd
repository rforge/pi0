\name{minque}
\alias{minque}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum norm quadratic unbiased estimation
}
\description{
This function computes the minimum norm quadratic unbiased estimate of variance components. Typically this is used as starting values of REML. 
}
\usage{
minque(y, varcov, start = rep(0, length(k)), lower.bound = -Inf, restricted = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A numeric vector of zero-mean response variable. 
}
  \item{varcov}{
	A list of variance-covariance matrices, the same as in \code{\link{varComp}} with no \code{random} argument. 
}
  \item{start}{
	A numeric vector of prior values of the ratio of variance components to the error variance. 
}
  \item{lower.bound}{
	A numeric value of the ratio of variance components to the error variance. If this is -Inf (default), minque solves the linear equation as in Rao (1972). If this is zero, quadratic programming is used to minimize the squared error between the two sides of the linear equation. 
}
  \item{restricted}{
	Not used. Currently only restricted likelihood is supported.
}
}

\value{A numeric vector of estimates of ratio of variance components to error variance. 
}
\references{
C. Radhakrishna Rao (1972) Estimation of Variance and Covariance Components in Linear Models. \emph{Journal of the American Statistical Association} 67, 112--115.
}
\author{
Long Qu
}

\seealso{
\code{\link{varComp.fit}}
}
\examples{
library(nlme)
data(Oxide)
vcf0 = varComp(Thickness~Source, Oxide, ~Lot/Wafer, 
      control=varComp.control(nlminb=nlminb.control(iter.max=0L)))
coef(vcf0,'var.ratio')
(st=minque(vcf0$residual.contrast, vcf0$working.cor, lower.bound = 0))
(vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer, control=varComp.control(start=st)))
coef(vcf,'var.ratio')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
