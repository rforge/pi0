\name{lassot.fit}
\alias{lassot.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting function for lassot regularized regression
}
\description{
The workhorse function for estimating the regression coefficients
}
\usage{
lassot.fit(x, y, lambdas = 10^seq(1, 5, length = 10), alphas = 1.5^seq(1, 17, length = 10), method = c("Coordinate", "Newton", "R"), eps = .Machine$double.eps^0.25 * ncol(x), niter = 1000L, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
}
  \item{alphas}{
%%     ~~Describe \code{alphas} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
