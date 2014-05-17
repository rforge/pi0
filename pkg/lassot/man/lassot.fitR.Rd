\name{lassot.fitR}
\alias{lassot.fitR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lassot.fitR(x.std, ym, lambdas, alphas, continuity = TRUE, eps = .Machine$double.eps^0.25 * ncol(x.std), niter = 1000L, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.std}{
%%     ~~Describe \code{x.std} here~~
}
  \item{ym}{
%%     ~~Describe \code{ym} here~~
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
}
  \item{alphas}{
%%     ~~Describe \code{alphas} here~~
}
  \item{continuity}{
%%     ~~Describe \code{continuity} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (x.std, ym, lambdas, alphas, continuity = TRUE, eps = .Machine$double.eps^0.25 * 
    ncol(x.std), niter = 1000L, verbose = TRUE) 
{
    if (isTRUE(continuity)) 
        alphas = pmax(alphas, continuityFactor * lambdas)
    nlambdas = length(lambdas)
    n = length(ym)
    p = ncol(x.std)
    betas = dfs = matrix(0, p, nlambdas)
    lastLambda = -Inf
    lastAlpha = -Inf
    for (i in seq(nlambdas)) {
        if (lambdas[i] == lastLambda) {
            betas[, i] = betas[, i - 1L]
            if (alphas[i] == lastAlpha) 
                next
        }
        else {
            lastLambda = lambdas[i]
            minAlphaDiff = suppressWarnings(min(abs(alphas[varComp::safeseq(i - 
                1L)] - alphas[i])))
            tmp.idx = which(abs(alphas[varComp::safeseq(i - 1L)] - 
                alphas[i]) == minAlphaDiff)
            tmp.i = which.min(abs(lambdas[tmp.idx] - lambdas[i]))
            if (length(tmp.i) > 0L) {
                betas[, i] = betas[, tmp.idx[tmp.i]]
            }
            else if (i > 1L) 
                betas[, i] = betas[, i - 1L]
        }
        lastAlpha = alphas[i]
        resids = ym - x.std \%*\% betas[, i]
        betas.old = betas[, i]
        for (iter in seq(niter)) {
            for (j in seq(p)) {
                if (TRUE) {
                  ans = lassot.fit1R(sum(x.std[, j] * resids) + 
                    betas[j, i], lambdas[i], alphas[i])
                  if (ans != betas[j, i]) {
                    resids = resids - x.std[, j] * (ans - betas[j, 
                      i])
                    betas[j, i] = ans
                  }
                }
                if (FALSE) {
                  if (betas[j, i] != 0) 
                    resids = resids + x.std[, j] * betas[j, i]
                  betas[j, i] = lassot.fit1R(sum(x.std[, j] * 
                    resids), lambdas[i], alphas[i])
                  if (betas[j, i] != 0) 
                    resids = resids - x.std[, j] * betas[j, i]
                }
                if (FALSE) {
                  beta.hat = xpy[j] - crossprod(xpx[j, ], betas[selected]) + 
                    betas[j, i]
                  ans = lassot.fit1R(beta.hat, lambdas[i], alphas[i])
                  if (ans != 0 && betas[j, i] == 0) {
                    selected = sort(c(selected, j))
                    xpx = crossprod(x.std, x.std[, selected, 
                      drop = FALSE])
                  }
                  betas[j, i] = ans
                }
            }
            mdiff = max(abs(betas.old - betas[, i]))
            if (isTRUE(verbose)) 
                cat("lambda=", lambdas[i], "\talpha=", alphas[i], 
                  "\titer=", iter, "\tmax.diff=", mdiff, "\n")
            if (mdiff <= eps) 
                break
            betas.old = betas[, i]
        }
        lasso.idx = (abs(betas[, i]) <= lassotFactor * alphas[i])
        dfs[lasso.idx, i] = as.numeric(betas[lasso.idx, i] != 
            0)
        for (j in which(!lasso.idx)) {
            bhat = sum(x.std[, j] * resids) + betas[j, i]
            A = -bhat
            B = alphas[i]^2 + lambdas[i] * alphas[i] * ffp1df
            C = -alphas[i]^2 * bhat
            tmpp = 2 * bhat^2 - 3 * alphas[i] * lambdas[i] * 
                ffp1df + 6 * alphas[i]^2
            tmpq = 8 * bhat^2 - (ffp1df * lambdas[i])^2 - 20 * 
                alphas[i] * lambdas[i] * ffp1df + 8 * alphas[i]^2
            M = 2 * A^3 - 9 * A * B + 27 * C
            K = A^2 - 3 * B
            N = M^2 - 4 * K^3
            sqrtN = sqrt(N)
            ninea2bqdsqrtN = 9 * alphas[i]^2 * bhat * tmpq/sqrtN
            dfs[j, i] = ((tmpp - ninea2bqdsqrtN)/cubeRoot(2 * 
                (M + sqrtN)^2) + (tmpp + ninea2bqdsqrtN)/cubeRoot(2 * 
                (M - sqrtN)^2) + 1)/3
        }
    }
    ans = betas
    attr(ans, "df") = colSums(dfs)
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
