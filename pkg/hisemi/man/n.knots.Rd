\name{n.knots}
\alias{directSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of spline knots
}
\description{
        Compute the number of spline knots to be the largest interger not greater than min(cutoff, max(0, n-cutoff)^rate).
}
\usage{
n.knots(n, cutoff=35, rate=0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ An integer vector of sample sizes. }
  \item{cutoff}{ A numeric vector of cutoff values. }
  \item{rate}{ A numeric vector of rates at which the number of knots increases with the sample size. }
}
\value{
 An integer vector, being the number of knots for splines. 
}
\references{ 

}
\author{ Long Qu \email{longorognol@hotmail.com}
}
%\note{ 
%}
%\seealso{ 
%}
\examples{
n.knots(10, 35, .2)
n.knots(35, 35, .2)
n.knots(135, 35, .2)
n.knots(1350, 35, .2)
n.knots(13500, 35, .2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
