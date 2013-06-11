\name{symbols}
\alias{plus}
\alias{colon}
\alias{star}
\alias{tilde}
%\docType{data}
\title{
	Some symbols that might be used in formulas
}
\description{
	Some symbols that might be used in formulas: "+", ":", "*", "~".
}

\usage{
plus
colon
star
tilde
}

\examples{
identical(plus, as.symbol("+"))
identical(colon, as.symbol(":"))
identical(star, as.symbol("*"))
identical(tilde, as.symbol("~"))
}
\keyword{datasets}
