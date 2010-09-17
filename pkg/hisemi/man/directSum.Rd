\name{directSum}
\alias{directSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Direct sum of matrices
}
\description{
        Direct sum of matrices, i.e., put matrices along the diagonal
}
\usage{
directSum(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ matrices to be \code{directSum}ed. }
}
\details{
The function treats each matrix as an element and put the element along the diagonal line. 
}
\value{
 A matrix, which is the direct sum of \dots. 
}
%\references{ 
%}
\author{ Long Qu \email{longor@iastate.edu}
}
%\note{ 
%}
%\seealso{ 
%}
\examples{
A=matrix(1:9,3,3)
all.equal(directSum(A), A) #TRUE
(B=directSum(-1,A))
directSum(A,B)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ math }
\keyword{ arith }
