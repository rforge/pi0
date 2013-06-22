\name{cholRoot}
\alias{cholRoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Lower Cholesky Root
}
\description{
	This function is similar to \code{\link[base:chol]{base::chol}}, but it returns the lower root and possibly remove unnecessary columns. 
}
\usage{
cholRoot(x)
}
\arguments{
  \item{x}{
	A positive semidefinite matrix. 
}
}
\details{
	This function calls \code{\link[base:chol]{base::chol}} with pivoting and returns the lower root. In the case of less than full rank \code{x}, unnecessary columns are removed. \code{tcrossprod(cholRoot(x))} will reproduce \code{x} (up to numerical errors). 
}
\value{
	A numeric matrix, not necessarily a square matrix. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link[Matrix:chol]{Matrix::chol}}, \code{\link[base:chol]{base::chol}}. 
}
\examples{

	set.seed(3456L)
	A=tcrossprod(matrix(runif(10L), 5L))
	cholRoot(A)
	range(A - tcrossprod(cholRoot(A)))
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ algebra }