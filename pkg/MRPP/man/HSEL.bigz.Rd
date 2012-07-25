\name{sample1}
\Rdversion{1.1}
\alias{sample1}
\alias{rejsample}
\alias{HSEL.bigz}
\title{
Sample without replacement when the population can be huge. 
}
\description{General a small sample of size M out of a huge population of size N. 
}
\usage{
HSEL.bigz(N, M, seed=0)
rejsample(N, M, seed=0)
sample1(N, seed=0)
}
\arguments{
  \item{N}{The size of the population.
}
  \item{M}{The sample size.
}
  \item{seed}{A \code{bigz} integer of random number seed. See \code{\link[gmp]{urand.bigz}}. 
}
}
\details{
\code{HSEL.bigz} generates a sample of size M from a population of size N, with space complexity linear in M. 
\code{rejsample} is similar to \code{HSEL.bigz}, but implemented using direct rejection through \code{\link{unique}}. 
\code{sample1} generates one number of 0 to N-1 using simple rejection and \code{\link[gmp]{urand.bigz}}. 
}
\note{
\code{rejsample} can be very slow when M is large. 
}
\value{
For \code{HSEL.bigz} and \code{rejsample}, a \code{bigz} vector of length \code{N}, between 1 and N. 
For \code{sample1}, a \code{bigz} integer, between 0 and N-1. 
}
\references{
Jarmo Ernvall and Olli Nevalainen. 1982. An Algorithm for Unbiased Random Sampling. THE COMPUTER JOURNAL, VOL. 25, NO. 1, 45--47
}
\seealso{
\code{\link{sample}}, \code{\link[gmp]{urand.bigz}}
}
\examples{
HSEL.bigz(factorialZ(15),5,12345)
rejsample(factorialZ(15),5,12345)
sample1(factorialZ(15), 12345)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{iteration}
\keyword{arith}