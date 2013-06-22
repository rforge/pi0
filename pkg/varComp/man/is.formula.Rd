\name{is.formula}
\alias{is.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Check if an object is a formula
}
\description{
	This function checks if the object inherits the 'formula' class. 
}
\usage{
is.formula(x)
}
\arguments{
  \item{x}{
	An R object.
}
}
\value{A logical scalar.
}
\author{
	Long Qu
}
\seealso{
	\code{\link[stats:formula]{stats::formula}}
}
\examples{
is.formula(~a+b+a:b)
is.formula('~a+b+a:b')

}
\keyword{ models}
