\name{factoradic}
\Rdversion{1.1}
\alias{dec2FR}
\alias{dec2permvec}
\alias{FR2dec}
\alias{FR2permvec}
\alias{permvec2dec}
\alias{permvec2FR}
\title{Factoradic numbers
}
\description{
Utility functions for factoradic numbers 
}
\usage{
FR2dec(FR)
dec2FR(dec,N)
FR2permvec(FR,base=1L)
permvec2FR(permvec)
dec2permvec(dec,N,base=1L)
permvec2dec(permvec)
}
\arguments{
  \item{FR}{Vector representation of factoradic number. 
}
  \item{dec}{Decimal representation of factoradic number.
}
  \item{N}{An integer}
  \item{base}{ 0 or 1, indicating the start index for a vector.}
  \item{permvec}{An integer vector of permutations of integers of size N.}
}
\details{
\code{FR2dec} and \code{dec2FR} convert factoradic numbers from/to decimal to/from vector representations. 
\code{FR2permvec} and \code{permvec2FR} converst factoradic numbers in vector representations from/to permutatin indices. 
\code{dec2permvec} and \code{permvec2dec} are convenient wrappers that converst factoradic numbers in vector representations from/to permutatin indices. 
}
\value{
For \code{FR2dec} and \code{permvec2dec}, an integer; 
for \code{dec2FR} and \code{permvec2FR}, a numeric vector;
for \code{FR2permvec} and \code{dec2permvec}, an integer vector.
}
\references{
Chapter 10 of .NET Test Automation Recipes

http://en.wikipedia.org/wiki/Factoradic

http://www.mathe2.uni-bayreuth.de/frib/KERBER/h00/node30.html
}
\seealso{
\code{\link{permuteTrt}}
}
\examples{
FR2dec(dec2FR(10,10))
set.seed(2340)
FR2permvec(permvec2FR(sample(10)))
permvec2dec(dec2permvec(100,10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{arith}
