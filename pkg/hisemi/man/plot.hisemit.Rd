\name{plot.hisemit}
\alias{plot.hisemit}
\alias{plotHisemitResid}
\alias{plotHisemitTuning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a hisemit object }
\description{
  Plot an object of class \code{hisemit}
}
\usage{
\method{plot}{hisemit}(x, type = c("tuning", "residual"), ...)
plotHisemitResid(obj, y.type = c("hist", "scatter"), 
        x.type = c("lfdr", "pi0", "f"), ...)
plotHisemitTuning(obj, SE = FALSE, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,obj}{ The object of class \code{hisemit} }
  \item{type}{ Either \code{'tuning'} or \code{'residual'}, the type of plot requested.  }
  \item{\dots}{ See \code{Details}.}

  \item{y.type}{ Either \code{'hist'} or \code{'scatter'}, the type of residual plot requested. }
  \item{x.type}{ One of \code{'lfdr'}, \code{'pi0'} or \code{'f'}, the x-axis when \code{y.type='scatter'}. }

  \item{SE}{ A logical scalar, indicating whether standard error bars should be added to the plot. }
  \item{add}{ A logical scaler, indicating whether a new plot should be generated or adding lines to the existing plot. }

}
\details{The generic function \code{plot.hisemit} calls either \code{plotHisemitResid} or \code{plotHisemitTuning} depending on \code{type}.
        For \code{plot.hisemit}, the \dots is the additinoal arguments to be passed to \code{plotHisemitResid} or \code{plotHisemitTuning}. 
        For residual plot (\code{plotHisemitResid}), the \dots is the additional arguments to be passed to \code{\link{hist}} when \code{y.type='hist'}; and to \code{\link{residuals.hisemit}} when \code{y.type='scatter'}.
        For tuning plot (\code{plotHisemitTuning}), the \dots is the additional arguments to be passed to either \code{\link{plot}} or \code{\link{lines}} depending on \code{add}. 
}
\value{
For histograms of residuals, an object from \code{\link{hist}} is returned. 
For scatter plot of residuals, an object from \code{\link{plot}} is returned. 
For tuning plot, a matrix with the range of tuning criterion is returned. 
}
\references{ Long Qu, Dan Nettleton, Jack Dekkers (2012) A hierarchical semiparametric model for incorporating inter-gene relationship information for analysis of genomic data. Biometrics, 68(4):1168-1177}
\author{ Long Qu \email{long.qu@wright.edu} }
%\note{ 
%}
\seealso{ \code{\link{plot.hisemit}}, \code{\link{fitted.hisemit}}, \code{\link{coef.hisemit}}, \code{\link{vcov.hisemit}}, \code{\link{residuals.hisemit}}, \code{\link{logLik.hisemit}}, \code{\link{confint.hisemit}}, \code{\link{print.hisemit}},
\code{\link{penLik.EMNewton}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}
