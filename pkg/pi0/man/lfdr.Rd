\name{lfdr and ppee}
\Rdversion{1.1}
\alias{lfdr}
\alias{lfdr.default}
\alias{lfdr.parncpt}
\alias{lfdr.sparncpt}
\alias{lfdr.nparncpt}
\alias{lfdr.nparncpp}
\alias{lfdr.CBUM}
\alias{lfdr.znormix}
\alias{lfdr.convest}
\alias{lfdr.discTMix}
\alias{ppee}
\alias{ppee.default}
\alias{ppee.parncpt}
\alias{ppee.sparncpt}
\alias{ppee.nparncpt}
\alias{ppee.nparncpp}
\alias{ppee.CBUM}
\alias{ppee.znormix}
\alias{ppee.convest}
\alias{ppee.discTMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Local False Discovery Rates (lfdr), i.e., Posterior Probability of Equivalent Expression (ppee) for gene expression profiling data. 
}
\description{
        These functions return pi0*(null test statistic density)/(marginal test statistic density). 
}
\usage{
lfdr(object, ...)
ppee(object, ...)
\method{lfdr}{default}(object, ...)
\method{ppee}{default}(object, ...)
\method{lfdr}{parncpt}(object, ...)
\method{ppee}{parncpt}(object, ...)
\method{lfdr}{sparncpt}(object, ...)
\method{ppee}{sparncpt}(object, ...)
\method{lfdr}{nparncpt}(object, ...)
\method{ppee}{nparncpt}(object, ...)
\method{lfdr}{nparncpp}(object, ...)
\method{ppee}{nparncpp}(object, ...)
\method{lfdr}{CBUM}(object, ...)
\method{ppee}{CBUM}(object, ...)
\method{lfdr}{znormix}(object, ...)
\method{ppee}{znormix}(object, ...)
\method{lfdr}{convest}(object, ...)
\method{ppee}{convest}(object, ...)
\method{lfdr}{discTMix}(object, ...)
\method{ppee}{discTMix}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of corresponding classes.
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{Other arguments currently not used.}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{ A numeric vector of lfdr. 
}
\author{Long Qu \email{lqu@biostatsolutions.com}
%%  ~~who you are~~
}
\note{
        \code{lfdr} and \code{ppee} are equivalent. 
%%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{parncpt}}, \code{\link{nparncpt}}, \code{\link{sparncpt}},  \code{\link{nparncpp}}, \code{\link{fdr}},
\code{\link{CBUM}}, \code{\link{znormix}}, \code{\link{convest}}, \code{\link{discTMix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
