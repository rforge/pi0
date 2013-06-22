\name{normalizeTrace}
\alias{normalizeTrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Rescales a square matrix such that the trace is the same as the dimension
}
\description{
	Rescales a square matrix such that the trace is the same as the dimension. 
}
\usage{
normalizeTrace(x)
}
\arguments{
  \item{x}{
	A square numeric matrix to be rescaled. 
}
}
\details{This is computed as \code{{x=as.matrix(x); x/mean(diag(x))}}. 
}
\value{
	A numeric matrix. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link[base:diag]{base::diag}}
}
\examples{
	set.seed(3456L)
	A=tcrossprod(runif(10L))
	A0 = normalizeTrace(A)
	sum(diag(A0)) - nrow(A)
	
}
\keyword{ array }
