\name{grid.search}
\Rdversion{1.1}
\alias{grid.search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs a grid search to minimize the objective function 
}
\description{Performs a grid search to minimize the objective function 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grid.search(obj, lower, upper, ngrid, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{objective function to be minimized
%%     ~~Describe \code{obj} here~~
}
  \item{lower}{numeric vector giving the lower bound of grid for each dimension
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{numeric vector giving the upper bound of grid for each dimension
%%     ~~Describe \code{upper} here~~
}
  \item{ngrid}{numeric vector giving the number of points each dimension
%%     ~~Describe \code{ngrid} here~~
}
  \item{\dots}{other arguments passed to \code{obj}
%%     ~~Describe \code{\dots} here~~
}
}
\details{This function first call \code{expand.grid} then evaluate \code{obj} to find a minimum. The number of calls to \code{obj} is \code{prod(ngrid)}. 
This is useful for finding a good starting values for many optimization routines. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{a numeric vector of the parameter that minimizes \code{obj}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{Long Qu \email{lqu@biostatsolutions.com}
%%  ~~who you are~~
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{optim}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function(obj, lower, upper, ngrid, ...)
%%{
%%    p=max(c(length(lower),length(upper),length(ngrid)))
%%    lower=rep(lower, length=p)
%%    upper=rep(upper, length=p)
%%    ngrid=rep(ngrid, length=p)
%%
%%    knot.list=list()
%%    for(i in 1:p) knot.list[[i]]=seq(from=lower[i], to=upper[i], length=ngrid[i])
%%    names(knot.list)=paste(names(lower),names(upper),names(ngrid),sep='.')
%%    
%%    grids=do.call(expand.grid, knot.list)
%%    ans=apply(grids, 1, obj, ...)
%%    return(unlist(grids[which.min(ans),]))
%%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
