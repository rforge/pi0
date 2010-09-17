\name{logit}
\alias{logit}
\alias{logit.inv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logit link and its inverse }
\description{
        Logit link and its inverse
}
\usage{
logit(mu)
logit.inv(eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{  See \code{\link{make.link}}. }
  \item{eta}{ See \code{\link{make.link}}. }
}
\details{
        \code{logit} is defined as \code{make.link("logit")$linkfun}. 
        \code{logit.inv} is defined as \code{make.link("logit")$linkinv}.
}
\value{
        A link function or its inverse. See \code{\link{make.link}}. 
}
%\references{ }
\author{ Long Qu \email{longor@iastate.edu}}
%\note{
%}
\seealso{\code{\link{make.link}}}
%\examples{
%}
\keyword{ models }
