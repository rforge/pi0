\name{Internal formula manipulation}
\alias{callList2terms}
\alias{splitTerm}
\alias{sortTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Internal formula term manipulation functions
}
\description{
	These are helper function used in the  formula interface of \code{\link{varComp}} and \code{\link{varComp.test}}. They either sort or split the formula terms involving ":". 
}
\usage{
callList2terms(cl)
splitTerm(term)
sortTerm(term, priority)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
		A list converted from a call. 
}
	\item{term}{
		A character scalar of a formula term to be re-ordered alphanumerically with respect to the components separated by ":".
}
	\item{priority}{
		A character vector whose elements have higher priority to be ordered to the front of the term. 
}
}
\details{
	The main use of \code{sortTerm} is to transform an interaction term like "B:A" into ordered form "A:B" such that redundant terms are easier to be identified. The results will depend on the system locale.
}
\value{
	\code{callList2terms} returns a character vector of individual terms (possibly with duplicates). \code{splitTerm} returns a character vector of unique individual components in the \code{term}. \code{sortTerm} returns a character scalar with individual components sorted. 
}
\author{
	Long Qu
}
\seealso{
	\code{\link[stats:formula]{stats::formula}}
}
\examples{
sortTerm("sex:ibs(SNP)")  ## ibs(SNP) is moved to the front in most locales
sortTerm("sex:ibs(SNP)", "sex")  ## the same as input
splitTerm("sex:ibs(SNP)")	## two components
sortTerm("ns(x, df=3):a:b:ibs(SNP)", "b")  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
