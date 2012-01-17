\name{plot.sparncpt}
\Rdversion{1.1}
\alias{plot.sparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot an object of class sparncpt, i.e., semiparametric estimate of noncentrality parameters
}
\description{Plot the histogram of observed t-statistics together with its fitted density estimate; also plotted is the estimated density of noncentrality parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{sparncpt}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{sparncpt}
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{ currently not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{ Left panel shows the density estimate of observed t-statistics, overlapped with a histogram; right panel shows the estimated density of noncentrality parameters. Solid line is the actual mean of the estimate; dashed line is located at zero. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{the \code{\link{invisible}} \code{x} itself
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
\examples{
\dontrun{
data(simulatedTstat)
(npfit=nparncpt(tstat=simulatedTstat, df=8)); 
(pfit=parncpt(tstat=simulatedTstat, df=8, zeromean=FALSE)); plot(pfit)
(pfit0=parncpt(tstat=simulatedTstat, df=8, zeromean=TRUE)); plot(pfit0)
(spfit=sparncpt(npfit,pfit)); plot(spfit)
}
}
% R documentation directory.
\keyword{ hplot }
\keyword{ methods }
