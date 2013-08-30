\name{scaledTMix.null}
\alias{scaledTMix.null}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit the null model to t-statistics }
\description{
  This function fits a null model to t-statistics, i.e., a two-component mixture, with one component being central t-distribution, the other component being scaled central t-distribution with scale parameter larger than 1. The mixing proportion for the central t-distribution is pi0. 
}
\usage{
scaledTMix.null(tstat, df, starts = 
        list(pi0 = seq(0.1, 0.99, length = 20), 
             scale = 2^seq(0.01, log2(max(abs(tstat))),
                length = 20)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A numeric vector of t-statistics. }
  \item{df}{ A numeric scalar or vector of the same length as \code{tstat}, giving the degrees of freedom for the \code{tstat}. }
  \item{starts}{ A list of two components, \code{pi0} and \code{scale}. Each being a numeric vector, which defines the grid for searching starting values. }
}
\details{
  The function finds maximum likelihood estimates of \code{pi0} and \code{scale}. \code{pi0} should lie between 0 and 1. \code{scale} should be larger than 1. The \code{L-BFGS-B} method is used in optimization function \code{\link{optim}}. 
  The \code{\link[pi0]{parncp}} function with \code{zeromean=TRUE} fits the same model, except that it returns \code{ncpest} object instead of \code{hisemit} object.
}
\value{
  A \code{hesimit} object. 
}
\references{
Qu, Nettleton, Dekkers (2012) A Hierarchical Semiparametric Model for Incorporating Inter-gene Information for Analysis of Genomic Data. Biometrics, 68(4):1168-1177
} 
\author{ Long Qu \email{long.qu@wright.edu}
}
%\note{ 
%  ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link[pi0]{parncp}}, \code{\link{scaledTMix.sat}}, \code{\link{scaledTMix.psat}}}
\examples{
set.seed(99927220)
(tstat=rt(5,1))
scaledTMix.null(tstat,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }% __ONLY ONE__ keyword per line
