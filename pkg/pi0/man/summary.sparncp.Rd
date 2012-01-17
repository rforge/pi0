\name{summary.sparncpt}
\alias{summary.sparncpt}
\alias{print.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{        Print summary for an object of class sparncp
 }
\description{
        Print summary for an object of class \code{sparncpt}. This includes the estimated proportion of true null hypotheses (\code{pi0}),
estimated mean of noncentrality parameters (\code{mu.ncp}), estimated standard deviation of noncentrality parameters (\code{sd.ncp}), and 
effective number of parameters (\code{enp}).
}
\usage{
\method{summary}{sparncpt}(object, ...)
\method{print}{sparncpt}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,x}{ an object of class \code{sparncpt}}
  \item{\dots}{ currently not used }
}
%\details{
%}
\value{ the \code{\link{invisible}} \code{object} itself
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Qu L, Nettleton D, Dekkers JCM. (2012) Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. Biometrics (in press).
}
\author{Long Qu 
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{sparncpt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
