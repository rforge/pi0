\name{bw.mse.pdf.asym}
\Rdversion{1.1}
\alias{bw.mse.pdf.asym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find asymptotically optimal bandwidth for KDE by minimizing MSE at the first data point
}
\description{
Find asymptotically optimal bandwidth for KDE by minimizing MSE at the first data point
}
\usage{
bw.mse.pdf.asym(x,iter=FALSE,eps=1e-3,iter.max=100,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{iter.max}{
%%     ~~Describe \code{iter.max} here~~
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x,iter=FALSE,eps=1e-3,iter.max=100,verbose=FALSE)
{require(ks)

        bw0=bw.nrd(x)
        Rkern=density(x,bw0,from=x[1],to=x[1],n=1,give.Rkern=TRUE)

        n.iter=1
        repeat{
            f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
            ddf=drvkde(x,2,bw0,se=FALSE)
            ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
            bw=(f.1/ddf.1/ddf.1*Rkern/length(x))^.2
            if(!iter || abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
            if(verbose)cat(bw,fill=T)
            bw0=bw
            n.iter=n.iter+1
        }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
