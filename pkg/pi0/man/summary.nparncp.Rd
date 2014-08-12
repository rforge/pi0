\name{summary.nparncpt}
\Rdversion{1.1}
\alias{summary.nparncpt}
\alias{print.nparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Print summary for an object of class nparncp
}
\description{
        Print summary for an object of class \code{nparncpt}. This includes the estimated proportion of true null hypotheses (\code{pi0}),
estimated mean of noncentrality parameters (\code{mu.ncp}), estimated standard deviation of noncentrality parameters (\code{sd.ncp}), 
the estimated effective number of parameters (\code{enp}), and the tuning parameter (\code{lambda}) that controls the amount of smoothing. 
}
\usage{
\method{summary}{nparncpt}(object, ...)
\method{print}{nparncpt}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,x}{ an object of class \code{nparncpt}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{ the \code{\link{invisible}} \code{object} itself
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Qu L, Nettleton D, Dekkers JCM. (2012) Improved Estimation of the Noncentrality Parameter Distribution from a Large Number of $t$-statistics, with Applications to False Discovery Rate Estimation in Microarray Data Analysis. Biometrics, 68, 1178--1187.
}
\author{ Long Qu 
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
