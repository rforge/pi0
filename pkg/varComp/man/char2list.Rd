\name{char2list}
\alias{char2list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a list with components given by their names
}
\description{
	For a given vector of object names, \code{char2list} creates a named list with all these objects as components. 
}
\usage{
char2list(Names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Names}{
	A character vector of object names.
}
}
\value{
	A named list of the same length as \code{Names}, each element being the corresponding object with that name. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link[base:get]{base::get}}
}
\examples{
char2list( c('letters', 'c') )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities	 }