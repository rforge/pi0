\name{   loadOrInstall       }
\alias{    loadOrInstall     }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Load packages or install and load packages
}
\description{
This functino loads packages if they are installed already. If failed, it tries to install the package automatically and re-load the package after installation. 
}
\usage{
loadOrInstall(package, dependencies="Depends",...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{package}{A character vector of the names of packages to be loaded. }
  \item{dependencies}{The same as in \code{\link{install.packages}}. }
  \item{\dots}{Additional arguments passed to \code{\link{require}}, \code{\link{install.packages}}. }
}
\details{
The function calls \code{\link{require}} to see if the package can be loaded. If not, it calls \code{\link{install.packages}}. 
Then \code{\link{require}} is called again for the installed package. 
}
\value{
A logical scalar. \code{TRUE} if all packages are loaded sucessfully. If \code{FALSE}, its \code{"failed"} attribute contains a vector of package names that failed to be loaded. 
}
%\references{
%
%
%}
\author{ Long Qu [\email{lqu@biostatsolutions.com}] }
%\note{
%
%
%}
\seealso{ 
    \code{\link{library}}, \code{\link{require}}, \code{\link{install.packages}}

}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities  }
%\keyword{  }
