\name{get.dp.dw.kde}
\Rdversion{1.1}
\alias{get.dp.dw.kde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute variable importance measure of partial derivative of p-values w.r.t. weights
}
\description{
Compute variable importance measure of partial derivative of p-values w.r.t. weights
}
\usage{
get.dp.dw.kde(y, perm.mat, r = 1:ncol(as.matrix(y)), test = FALSE, distObj = dist(y), 
	mrpp.stats = mrpp.test.dist(distObj, perm.mat = perm.mat,
	wtmethod=wtmethod[1])$all.stat, bw = bw.mse.pdf.asym(mrpp.stats), 
	cperm.mat, wtmethod=integer(1), min.wts = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{perm.mat}{
%%     ~~Describe \code{perm.mat} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
}
  \item{distObj}{
%%     ~~Describe \code{distObj} here~~
}
  \item{mrpp.stats}{
%%     ~~Describe \code{mrpp.stats} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{cperm.mat}{
%%     ~~Describe \code{cperm.mat} here~~
}
  \item{wtmethod}{
%%     ~~Describe \code{cperm.mat} here~~
}
  \item{min.wts}{
%%     ~~Describe \code{min.wts} here~~
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
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
