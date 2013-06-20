\name{internals}
\alias{AEI}
\alias{AOI}
\alias{EI}
\alias{OI}
\alias{WAI}
\alias{PREML}
\alias{preprocPREML}
\alias{obj}
\alias{obj2}
\alias{score}
\alias{updateLI}
\alias{updateLIkLI}
\alias{updateLIy}
\alias{updateNegGrad}
\alias{updateNegHess}
\alias{updateNums}
\alias{updateNums2}
\alias{updateNumsPart}
\alias{updateTr1}
\alias{updateTr2}
\alias{V}
\alias{updateDenom}
\alias{hess}
\alias{gradFunc}
\alias{varComp.test.2modelDoTest}
\alias{varComp.test.nulDoTest}
\alias{varComp.test.altDoTest}
\alias{varComp.test.Common}
\alias{varComp.LinScore.Normal}
\alias{varComp.LinScore.Satterthwaite}
\alias{varComp.LinScore.SSAS155}
\alias{varComp.LinScore.test}
\alias{varComp.RLRT.test}
\alias{varComp.VM03.test}
\alias{varComp.SS95.test}

\title{	Internal functions not supposed to be called directly
}
\description{
These internal functions are mainly used by \code{varComp.fit} and \code{varComp.test.Common}. These functions make reference to and make changes to non-local objects. Thus, they are not supposed to be called directly, unless their environments are set properly.
}
\usage{
AEI()
AOI()
EI()
OI()
WAI()
PREML()
preprocPREML(tau)
obj(tau)
obj2(ltau)
score()
updateLI(tau)
updateLIkLI()
updateLIy()
updateNegGrad()
updateNegHess()
updateNums()
updateNums2()
updateNumsPart()
updateTr1()
updateTr2()
V(tau)
updateDenom()
hess(tau)
gradFunc(tau)

varComp.test.nulDoTest(null.fit, additional.varcov, test = "LinScore", 
                       control = varCompTest.control(test), alt.fit = NULL, ...)
varComp.test.2modelDoTest(null.fit, alt.fit, test = "LinScore", 
                          control = varCompTest.control(test), ...)
varComp.test.altDoTest(alt.fit, null = integer(0L), test = "LinScore", 
                       control = varCompTest.control(test), null.fit = NULL, ...)
varComp.test.Common()


varComp.LinScore.Normal(all.scores, lin.form, infoMat, null, w, tr1, n, LIkLI, 
                        tau.idx, ...)
varComp.LinScore.Satterthwaite(all.scores, infoMat, null, w, tr1, n, ...)
varComp.LinScore.SSAS155(all.scores, lin.form, infoMat, null, w, tr1, n, 
                           LIkLI, tau.idx, non.pd, control, ...)

varComp.LinScore.test(control, n, tau.idx, LIkLI, tr1, infoMat, all.scores, non.pd)
varComp.RLRT.test(control, negHess, alt.fit, null.fit, tau.idx)
varComp.VM03.test(control, infoMat, tau.idx, LIkLI, tr1, n, LIy, all.scores)
varComp.SS95.test(control, infoMat, tau.idx, LIkLI, LIy, tr1, n, all.scores)

}
\arguments{
\item{
	tau,ltau,null.fit,additional.varcov,test,control,
	alt.fit,...,null,all.scores,lin.form,infoMat,w,tr1,
	n,LIkLI,tau.idx,non.pd,negHess,LIy
}{To be described later, if indeed needed.}
}
%- maybe also 'usage' for other objects documented here.
% \details{
%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
% ...
% }
% \references{
% ~put references to the literature/web site here ~
% }
% \author{
%  ~~who you are~~
% }
% \note{
%  ~~further notes~~
% }

% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function () 
% {
    % 0.5/(n + 2) * (n * n * nums2/denom/2 - outer(tr1, tr1))
  % }
% }
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
