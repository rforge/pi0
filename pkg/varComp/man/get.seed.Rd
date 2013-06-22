\name{get.seed}
\alias{get.seed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Recording pseudo-random number seeds
}
\description{
\code{get.seed} obtains current pseudo-random number seeds. 
}
\usage{
get.seed()
}
%- maybe also 'usage' for other objects documented here.
\details{
This function obtain the \code{.Random.seed} object in the global environment. If it is absent, \code{runif(1L)} is called and then the seed is obtained. 
}
\value{
	A numeric vector of \code{.Random.seed}, with \code{RNGkind} attribute being the result from calling \code{RNGkind()}.
}
\author{
	Long Qu
}
\seealso{
	\code{\link[base:Random]{base::.Random.seed}}
}
\examples{
	set.seed(2034L)
	all.equal(get.seed(), .Random.seed, check.attributes = FALSE)
}
\keyword{ distribution }
