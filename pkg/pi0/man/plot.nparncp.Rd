\name{plot.nparncpt}
\Rdversion{1.1}
\alias{plot.nparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot an object of class nparncpt, i.e., nonparametric estiamte of noncentrality parameters
}
\description{Plot the Network information criterion (NIC), effective number of parameters (ENP), and estimated proportion (pi0) of true null hypotheses
for different choices of tuning parameters; also plot the estimated density of noncentrality parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{nparncpt}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{nparncpt}
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{currently not used. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{For NIC, only values within 2 s.e.'s of the minimum are shown. The solid line on NIC, ENP and pi0 shows the final tuning parameter, i.e., the one that minimizes NIC. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{Invisible \code{par}.
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

\seealso{\code{\link{nparncpt}}, \code{\link{sparncpt}},\code{\link{parncpt}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
data(simulatedTstat)
(npfit=nparncpt(tstat=simulatedTstat, df=8, plotit=FALSE)); plot(npfit)
(pfit=parncpt(tstat=simulatedTstat, df=8, zeromean=FALSE)); plot(pfit)
(pfit0=parncpt(tstat=simulatedTstat, df=8, zeromean=TRUE)); plot(pfit0)
(spfit=sparncpt(npfit,pfit)); plot(spfit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot}
\keyword{ methods }
