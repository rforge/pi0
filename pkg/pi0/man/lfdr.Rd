\name{lfdr and ppee}
\Rdversion{1.1}
\alias{lfdr}
\alias{lfdr.parncpt}
\alias{lfdr.sparncpt}
\alias{lfdr.nparncpt}
\alias{ppee}
\alias{ppee.parncpt}
\alias{ppee.sparncpt}
\alias{ppee.nparncpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Local False Discovery Rates (lfdr), i.e., Posterior Probability of Equivalent Expression (ppee) for gene expression profiling data. 
}
\description{
        These functions return pi0*(central t-density)/(marginal t-density). 
}
\usage{
lfdr(object)
ppee(object)
\method{lfdr}{parncpt}(object)
\method{ppee}{parncpt}(object)
\method{lfdr}{sparncpt}(object)
\method{ppee}{sparncpt}(object)
\method{lfdr}{nparncpt}(object)
\method{ppee}{nparncpt}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{parncpt}, \code{sparncpt}, or \code{nparncpt}.
%%     ~~Describe \code{object} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{ A numeric vector of lfdr. 
}
\author{Long Qu \email{longor@iastate.edu}
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{sparncpt}}, \code{\link{fdr}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
