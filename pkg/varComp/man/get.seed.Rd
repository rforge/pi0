\name{get.seed}
\alias{get.seed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Pseudo-random number seeds
}
\description{
\code{get.seed} obtains current pseudo-random number seeds. 
}
\usage{
get.seed()
}
%- maybe also 'usage' for other objects documented here.
\details{
This function obtain the '.Random.seed' object in the global environment. If it is absent, \code{runif(1L)} is called and then the seed is obtained. 
}
\value{
	A numeric vector of '.Random.seed'. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link{.Random.seed}}
}
\examples{
	set.seed(2034L)
	identical(get.seed(), .Random.seed)
}
\keyword{ sysdata }
