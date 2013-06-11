\name{safeseq}
\alias{safeseq}
\title{
	A ``safer'' seq
}
\description{
	This is the same as \code{seq(from , to, by)} except that when the sign of by is wrong, \code{safeseq} returns a zero-length vector rather than throwing an error. 
}
\usage{
safeseq(from = 1L, to = 1L, by = 1L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{The same as the \code{from} in \code{\link[base]{seq}}. 
}
  \item{to}{The same as the \code{to} in \code{\link[base]{seq}}. 
}
  \item{by}{The same as the \code{by} in \code{\link[base]{seq}}. 
}
  \item{\dots}{The same as in \code{\link[base]{seq}}. 
}
}
\value{
A vector of the same class as \code{by*(from-to)}, but with length being possibly zero. 
}
\author{
	Long Qu
}

\seealso{
\code{\link[base]{seq}}
}
\examples{
identical(integer(0L), safeseq(1L, 0L, 1L))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
