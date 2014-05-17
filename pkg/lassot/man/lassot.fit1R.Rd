\name{lassot.fit1R}
\alias{lassot.fit1R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lassot.fit1R(beta.hat, lambda, alpha, df = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.hat}{
%%     ~~Describe \code{beta.hat} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (beta.hat, lambda, alpha, df = FALSE) 
{
    pls = function(beta) 0.5 * (beta - beta.hat)^2 + lambda * 
        ifelse(abs(beta) <= lassotFactor * alpha, abs(beta), 
            ffp1df/2 * alpha * log(1 + beta^2/alpha^2))
    root0 = sign(beta.hat) * pmax(abs(beta.hat) - lambda, 0)
    if (abs(beta.hat) <= lambda + lassotFactor * alpha) {
        if (isTRUE(df)) 
            attr(root0, "df") = as.numeric(root0 != 0)
        return(root0)
    }
    A = -beta.hat
    B = alpha^2 + lambda * alpha * ffp1df
    C = -alpha^2 * beta.hat
    if (FALSE) {
        sols = polyroot(c(C, B, A, 1))
        sol = Re(sols[abs(Im(sols)) < 1e-06])
        obj = pls(sol)
        ans = sol[which.min(obj)]
        return(ans)
    }
    K = A * A - 3 * B
    M = 2 * A^3 - 9 * A * B + 27 * C
    N = M * M - 4 * K^3
    sqrtN = sqrt(N)
    if (N > 0) {
        ans = -oneThird * (A + cubRoot((M + sqrtN)/2) + cubRoot((M - 
            sqrtN)/2))
        if (isTRUE(df)) {
            tmpp = 2 * beta.hat^2 - 3 * alpha * lambda * ffp1df + 
                6 * alpha^2
            tmpq = 8 * beta.hat^2 - (ffp1df * lambda)^2 - 20 * 
                alpha * lambda * ffp1df + 8 * alpha^2
            ninea2bqdsqrtN = 9 * alpha^2 * beta.hat * tmpq/sqrtN
            attr(ans, "df") = ((tmpp - ninea2bqdsqrtN)/cubeRoot(2 * 
                (M + sqrtN)^2) + (tmpp + ninea2bqdsqrtN)/cubeRoot(2 * 
                (M - sqrtN)^2) + 1)/3
        }
        return(ans)
    }
    browser()
    theta = acos(-M/2/K/sqrt(K))
    ell = 0:2
    insideTrig = (theta + 2 * ell * pi)/3
    sol = 2/3 * sqrt(K) * cos(insideTrig) + oneThird * beta.hat
    if (FALSE) {
        dKdBetahat = 2 * beta.hat
        dThetadBetahat = -3 * alpha * ffp1df * beta.hat/4/K^2.5/sqrt(1 - 
            (alpha * ffp1df * lambda)^2/16/K^3)
        dSoldBetahat = (3 * cos(insideTrig) * dKdBetahat - 2 * 
            K * dThetadBetahat * sin(insideTrig) + 3 * sqrt(K))/9/sqrt(K)
    }
    obj = pls(sol)
    ans = sol[which.min(obj)]
    attr(ans, "obj") = obj
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
