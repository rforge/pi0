\name{subt}
\alias{subt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Subsampling a Microarray Data Set for Estimating Proportion of True Null Hypotheses }
\description{
This function subsamples the columns (arrays) of a microarray data set and do two-sample t-tests. Subsamples 
from each treatment group are obtained and combined. A t-test is conducted for each row (gene) of the 
subsampled data set and the p-value density at one is estimated for each combined subsample. 
}
\usage{
subt(dat, n1 = round(ncol(dat)/2), n2 = ncol(dat) - n1, 
      f1method = c("lastbin", "qvalue"), 
        max.reps = if(balanced)20 else 5, balanced = FALSE,  ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a numeric matrix, the microarray data set with each row being a gene, and each column being a 
                subject. The first \code{n1} columns correspond to treatment group 1 and the rest \code{n2} 
                columns correspond to treatment group 2. }
  \item{n1}{ a positive integer, the original sample size in treatment group 1. }
  \item{n2}{ a positive integer, the original sample size in treatment group 2. }
  \item{f1method}{ character, the name of the function to be used to estimate the p-value density at 1. The first argument of the function needs to be a vector of values. }
  \item{max.reps}{ a positive integer, the maximum number of subsamples to obtain per subsample size 
                configuration. If this is set to \code{Inf}, then all possible subsamples will be tried. 
                However, see Notes and the \code{R} argument of \code{\link{combn2R}}. }
  \item{balanced}{ logical, indicating whether only balanced subsamples are obtained. This is computationally 
                faster and is good for initial exploration purposes. }
%  \item{totalOnly}{logical. If \code{TRUE}, then only the total number of subsamples are returned, without
%                doing any actual subsampling. This is used only to see if the computational burden look 
%                reasonable or not before doing any actual computations. }
  \item{\dots}{ additional arguments used by \code{f1method}. }
}
\details{
This function tries to get possible subsamples through \code{\link{combn2R}}. \cr
For each total subsample size M=3,4,...,N, where N=n1+n2, do the following, 
\itemize{
\item{1}{For each treatment 1 subsample size m1=1,2,...,n1, let m2=M-m1. If 1<=m2<=n2 and at least one of \code{balanced} and m1=m2 is true, then do the following, 
    \itemize{
        \item{1.1}{Randomly choose \code{max.reps} subsamples among all possible subsamples by choosing m1 subjects from treatment group 1 and m2 subjects from treatment group 2, by using the function \code{\link{combn2R}} with \code{sample.method="diff2"} and \code{try.rest=TURE}. Note that this may \emph{not} be always possible due to some pratical computational limitations. See \code{\link{combn2R}} for details.}
        \item{1.2}{For each subsample obtained in \code{1.1},  (1) do a t-test for each gene (i.e., each row of the subsample), and (2) estimate the p-value density at one.}
    }
    }
}
}
\value{
%If \code{totalOnly=TRUE}, only a positive integer is returned, giving the total number of subsamples;\cr
%otherwise, 
an object of class \code{c("subt","matrix")}, which is a G-by-3 numeric matrix, where G is \code{nrow{dat}}, 
with column names 'f1', 'n1', and 'n2', corresponding to the p-value density at 1 and subsample size 
in each treatment group. This object also has the following \code{\link{attributes}},
\item{n1}{the same as the argument \code{n1}.}
\item{n2}{the same as the argument \code{n2}.}
\item{f1method}{the same as the argument \code{f1method}.}
\item{max.reps}{the same as the argument \code{max.reps}.}
\item{balanced}{the same as the argument \code{balanced}.}

}
\references{
Qu, L., Nettleton, D., Dekkers, J.C.M. Subsampling Based Bias Reduction in Estimating the Proportion of 
Differentially Expressed Genes from Microarray Data. Unpublished manuscript.
}
\author{ Long Qu \email{longor@iastate.edu} }
\note{ 
\code{max.reps} applies to each subsample size configuration. For example, 2 subjects subsampled from 
treatment group1 and 3 subjects subsampled from treatment group 2 will be considered as a different 
subsample size configuration than 3 subjects subsampled from treatment group 1 and 2 subjects subsampled 
from treatment group 2. For the small sample sizes commonly seen in microarray data, a large 
\code{max.reps} is rarely a big computational burden. But be careful when you do have a very large 
sample size, as the number of all possible subsamples grows very fast.
}
\seealso{\code{\link{print.subt}}, \code{\link{plot.subt}}, \code{\link{extrp.pi0}}, 
\code{\link{matrix.t.test}},\code{\link{combn2R}}, \code{\link{subex}}, \code{\link{lastbin}}, 
\code{\link[qvalue]{qvalue}}
}

\examples{
\dontrun{
set.seed(9992722)
## this is how the 'simulatedDat' data set in this package generated
simulatedDat=sim.dat(G=5000)        
## this is how the 'simulatedSubt' object in this package generated
simulatedSubt=subt(simulatedDat,balanced=FALSE,max.reps=Inf) 
}
data(simulatedSubt)
print(simulatedSubt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{iteration}
