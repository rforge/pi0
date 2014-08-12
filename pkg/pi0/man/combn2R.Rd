\name{combn2R}
\alias{combn2R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Randomly Choosing R Combinations of Two Groups of n and n2 Elements Taken m and m2 at a Time }
\description{
This is the enhanced version of \code{\link[utils]{combn}}. There are two groups, 
with n1 and n2 elements, respectively. Each time, m1 elements will be randomly chosen 
from group 1; and m2 elements will be randomly chosen from group 2. These m1+m2 elements 
form one combination. This function generate either all \code{choose(n1,m1)*choose(n2,m2)} 
such combinations or a subset of \code{R} of them. A function, possibly identity, is then 
applied to each selected combination. 
}
\usage{
combn2R(x, m, x2, m2, R, FUN = c, simplify = TRUE, 
        sample.method="all", try.rest=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ group 1 vector for combinations, or positive integer \eqn{n}{n} for \code{x <- seq(n)}. }
  \item{m}{ number of elements to choose from \code{x}, i.e., group 1. }
  \item{x2}{ group 2 vector for combinations, or positive integer \eqn{n_2}{n2} for \code{x2 <- seq(n)}. 
	     If missing, it reduces to get R combinations from x taken m at a time. See details. }
  \item{m2}{ number of elements to choose from \code{x2}, i.e., group 2. 
	     If missing, it reduces to get R combinations from x taken m at a time. See details. }
  \item{R}{ the number of combinations to be randomly chosen from all \code{choose(n,m)*choose(n2,m2)} 
	combinations. If missing or larger than all possible combinations, results from all combinations 
	will be \emph{tried}, but \emph{not} guaranteed when the total number of possible combinations 
	is too large. See details. }
  \item{FUN}{ a function to be applied to each chosen combination. When neither \code{x2} nor \code{m2} is missing, 
	and neither \code{m=0} nor \code{m2=0}, this function needs to accept at least two arguments, the first 
	of which is a vector of length \code{m}, which is a subset of \code{x}; and the second argument of 
	\code{FUN} is a vector of length \code{m2}, which is a subset of \code{x2}. Additional arguments are 
	supplied with \code{dots}. 
	When \code{combn2R} is only used for one group situtation (similar to \code{\link[utils]{combn}}),
	the second argument of \code{FUN} is not required. }
  \item{simplify}{ logical, indicating if the result should be simplified to an array (typically a matrix); 
	see \code{\link[utils]{combn}}. }
  \item{sample.method}{ character, specifying how samples should be taken, if not all combinations are 
	  generated; possible choices are "diff2", "all", and "replace" for two-group situation, and "replace" 
	  and "noReplace" in one-group situation. See Details. }
  \item{try.rest}{ logical, together with \code{sample.method}, specifying how samples should be taken, 
	if not all combinations are generated; see Details.  }
  \item{\dots}{ optionally, further arguments to \code{FUN}. }
}
\details{
 This function enhances \code{\link[utils]{combn}} in two ways.  
One is to deal with two-group situation, which is commonly seen in real designs; 
the other is to choose only a random sample of size \code{R} from all possible 
combinations to avoid unnecessary computation. \cr
When neither \code{x2} nor \code{m2} is missing and neither \code{m=0} nor \code{m2=0}, 
the function works in two-group mode. In this situation, \cr

(A)  if \code{sample.method="diff2"},  \code{combn2R} will try to sample \code{R} combinations from each 
group separately. That is, first sample \code{R} combinations from \code{x} taken \code{m} at a time, 
and then sample \code{R} combinations from \code{x2} taken \code{m2} at time. The results are then combined 
to give \code{R} combinations from the two groups. This sampling method will make the samples as different 
from each other as possible. But when \code{R} is larger than \code{min(choose(n1,m1),choose(n2,m2))}, 
it is not possible to get \code{R} samples from each group separately. If this happens and \code{try.rest=FALSE}, 
then \code{R} will be reset to \code{min(choose(n,m),choose(n2,m2))} and the function works as before;
otherwise, if \code{try.rest=TRUE}, then \code{sample.method} will be reset to "all" and the function will 
try to get \code{R} samples from all \code{choose(n1,m1)*choose(n2,m2)} combinations (see below). \cr

(B)  if \code{sample.method="all"},  \code{combn2R} will try to sample \code{R} samples from all
\code{choose(n1,m1)*choose(n2,m2)} combinations directly. This means two samples of size \code{m+m2} 
may have the same sample of size \code{m} (or \code{m2}) which comes from \code{x} (or \code{x2}). For 
example, if \code{x=1:3}, \code{m=1}, \code{x2=4:6}, \code{m=2} and \code{R=2}, then it is possible 
to get one sample to be \code{1} and \code{4,5}, but the other sample is \code{1} and \code{5,6}. That is, 
the same sample from \code{x} is used in both results. This will not happen when \code{sample.method='diff2'}.
However, this will guarantee any two samples of size \code{m+m2} will differ in at least one element. \cr

(C)  if \code{sample.method='replace'},  \code{combn2R} will not guarantee the uniqueness of the 
\code{R} combinations in any way. It is possible to have two exactly the same samples of size 
\code{m+m2}. \cr


Because the number of possible combinations grows very fast, computational limitations may be reached. In
this case, if \code{try.rest=TRUE}, then \code{sample.method} will be reset to "replace", which uses the 
least computational resources; otherwise, an error will be generated. \cr 

When either \code{x2} or \code{m2} is missing, or one of \code{m} and \code{m2} is zero, 
the function works in one-group mode. In this situation, \code{sample.method="diff2"} and 
\code{sample.method="all"} will be treated the same as \code{sample.method="noReplace"}, and \code{combn2R}
will try to obtain \code{R} different combinations from all possible combinations for the non-missing group. 
Again, if this fails due to computational limitations, \code{sample.method} will be reset to "replace" and 
no guarantee will be made to ensure the \code{R} combinations are different from each other. \cr



}

\value{
a \code{\link[base]{list}} or \code{\link[base]{array}} (in nondegenerate cases), similar to \code{\link[utils]{combn}}.
An attribute \code{"sample.method"} will be added to the \code{list} or \code{array}, which stores the 
\emph{actual} sampling method that has been used, which may or may not be the same as specified in the 
argument. 
}

\note{
Note that the results are \emph{not} necessarily in order, which is a difference from \code{\link{combn}}. 
}


\author{ Long Qu   }
\seealso{ \code{\link[utils]{combn}} in \pkg{utils} or \pkg{combinat},  }
\examples{
combn2R(4,2) ### the same as combn(4,2), except an additional attribute
combn2R(1:4*2,2)
combn2R(4,2,5,1)
combn2R(4,2,5,1,FUN=sum)
set.seed(992722)
combn2R(4,2,R=3) ### the same as combnR(4,2,3), except an additional attribute
combn2R(4,2,R=10) ### only 6 combinations possible
combn2R(4,2,5,1,R=8)
combn2R(1:4*2,2,5,1,R=50) ### only 30 combinations possible
combn2R(1:4*2,2,5,1,R=5) ### when considering only one group, there are several common samples.
### no common samples, even considering only one group
combn2R(1:4*2,2,5,1,R=5, sample.method="diff2") 
combn2R(1:3*3,1,3,1,R=5, sample.method="replace") ### two pairs of exactly common samples
combn2R(100,3,100,3,R=5, sample.method="all") ### 'all' combinations not feasible (~3e10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ iteration }
% __ONLY ONE__ keyword per line
