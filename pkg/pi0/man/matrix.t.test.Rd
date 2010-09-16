\name{matrix.t.test}
\alias{matrix.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply a Two Sample T-test to Each Row or Column of a Matrix }
\description{
This function applies the two sample t-test to each row or column of a matrix. 
}
\usage{
matrix.t.test(x, MARGIN = 1, n1 = if (MARGIN == 1) floor(ncol(x)/2)
    else floor(nrow(x)/2), n2 = if (MARGIN == 1) ncol(x) - n1 else 
    nrow(x) - n1, pool = TRUE, pOnly=TRUE, tOnly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric \code{matrix} to which the t-test will be applied to, by row or by column. }
  \item{MARGIN}{ either 1 or 2. If \code{MARGIN=1}, apply the t-test to each row of \code{x}; otherwise, 
		if \code{MARGIN=2}, apply the t-test to each column of \code{x}. See also \code{\link{apply}}.}
  \item{n1}{ sample size of the first group. It should be smaller than the appropriate \code{dim} of \code{x}.}
  \item{n2}{ sample size of the second group.  It should be smaller than the appropriate \code{dim} of \code{x}.}
  \item{pool}{ logical, indicating if the variance estimate should be pooled. If \code{FALSE}, 
	Welch (i.e. Satterthwaite) approximation to the degrees of freedom is used. }
  \item{pOnly}{ logical, indicating if only a vector of p-values should be returned. }
  \item{tOnly}{ logical, indicating if only a vector of t-statistics should be returned. 
		This argument overwrites \code{pOnly}, if they are conflicting. }
}
\details{
	This is a much faster function for "almost" the same purpose of \code{apply} each \code{MARGIN} of 
	\code{x} a \code{t.test}, i.e., the mean of the first \code{n1} elements is compared with the 
	mean of the rest \code{n2} elements, for each row or column depending on the \code{MARGIN}.
	See the Value section for differences. 

}
\value{
If \code{pOnly=TRUE} (the default situation), a numeric vector of p-values is returned, the length of which is determined by \code{MARGIN}. \cr
If \code{tOnly=TRUE} , a numeric vector of t-statistics is returned, the length of which is determined by \code{MARGIN}. \cr
If \code{tOnly=TRUE} \emph{and} \code{tOnly=TRUE}, a numeric vector of t-statistics is returned, the length of which is determined by \code{MARGIN}, as \code{tOnly} overwrites \code{pOnly}. \cr
If \code{pOnly=FALSE} \emph{and} \code{tOnly=FALSE}, a list of three components is returned:\cr
  \item{stat}{a numeric vector of the t-statistics, one for each row or column, depending on \code{MARGIN}.}
  \item{df}{a numeric vector of degrees of freedom. If \code{pool} is \code{TRUE}, this vector is of length 1, i.e. n1+n2-2; if \code{pool} is \code{FALSE}, this vector is of the same length as \code{stat}, depending on \code{MARGIN}.}
  \item{p.value}{a numeric vector of p-values, one for each row or column, depending on \code{MARGIN}.}
}

\author{ Long Qu \email{longor@iastate.edu}}
\seealso{\code{\link{apply}},\code{\link{t.test}} }
\examples{
set.seed(9992722)
dat=matrix(rnorm(30),3,10)
(pvals=matrix.t.test(dat,1,5,5)) # [1] 0.2112825 0.8366920 0.2891014
(pvals2=apply(dat,1,function(xx)t.test(xx[1:5],xx[6:10],var.equal=TRUE)$p.val))
all.equal(pvals,pvals2) ## TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
