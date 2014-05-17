\name{lassot.fit}
\alias{lassot.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lassot.fit(x, y, lambdas = 10^seq(1, 5, length = 10), alphas = 1.5^seq(1, 17, length = 10), method = c("Coordinate", "Newton", "R"), eps = .Machine$double.eps^0.25 * ncol(x), niter = 1000L, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
}
  \item{alphas}{
%%     ~~Describe \code{alphas} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (x, y, lambdas = 10^seq(1, 5, length = 10), alphas = 1.5^seq(1, 
    17, length = 10), method = c("Coordinate", "Newton", "R"), 
    eps = .Machine$double.eps^0.25 * ncol(x), niter = 1000L, 
    verbose = TRUE) 
{
    continuity = TRUE
    x = as.matrix(x)
    p = ncol(x)
    n = length(y)
    ym = as.double(y - mean(y))/sd(y)
    xm = sweep(x, 2L, colMeans(x), "-")
    xnorm = sqrt(colSums(xm * xm))
    x.std = sweep(xm, 2L, xnorm, "/")
    beta0 = drop(crossprod(x.std, ym))
    beta0 = beta0 * drop(crossprod(beta0)/crossprod(x.std \%*\% 
        beta0))
    grids = expand.grid(rev(sort(alphas)), rev(sort(lambdas)))
    lambdas = grids[, 2]
    alphas = grids[, 1]
    if (isTRUE(continuity)) 
        alphas = pmax(continuityFactor * lambdas, alphas)
    lambda.n = length(lambdas)
    method = match.arg(method)
    if (method == "Coordinate") {
        ans = .C("lassot", as.double(x.std), ym, length(ym), 
            p, b = rep(0, length = p * lambda.n), as.double(lambdas), 
            as.double(alphas), length(lambdas), as.double(eps), 
            as.integer(niter), as.integer(verbose))
    }
    else if (method == "R") {
        ans = lassot.fitR(x.std, ym, lambdas, alphas, continuity, 
            eps, niter, verbose)
        ans = list(b = ans)
    }
    else if (method == "Newton") {
        xpx = crossprod(x.std)
        obj = function(betas, lambda, alpha) 0.5 * sum(crossprod(ym - 
            x.std \%*\% betas)) + lambda * sum(ifelse(abs(betas) <= 
            alpha, log(1 + betas * betas/alpha/alpha), abs(betas)/a - 
            1 + log(2)))
        der = function(betas, lambda, alpha) -beta0 + xpx \%*\% 
            betas + lambda * ifelse(abs(betas) <= alpha, 2 * 
            betas/(betas * betas + alpha * alpha), ifelse(betas > 
            0, 1, -1)/alpha)
        ret = matrix(mean(y), p + 1, lambda.n)
        for (i in 1:lambda.n) ret[-1, i] = nlminb(rep(0, p), 
            obj, der, lambda = lambdas[i], alpha = alphas[i], 
            control = list(trace = 1))$par
    }
    betas = matrix(ans$b, p, lambda.n)
    sse = colSums((ym - x.std \%*\% betas)^2)
    gcv = sse/n/(1 - attr(ans$b, "df")/n)^2
    ibest = which.min(gcv)
    ret = c(mean(y), betas[, ibest]/xnorm * sd(y))
    attr(ret, "alpha") = alphas[ibest]
    attr(ret, "lambda") = lambdas[ibest]
    attr(ret, "sig2e") = crossprod(y - x \%*\% ret[-1] - mean(y))/(n - 
        min(n - 1, attr(ans$b, "df")[ibest]))
    attr(ret, "df") = attr(ans$b, "df")[ibest]
    attr(ret, "gcv") = gcv
    attr(ret, "betas") = betas/xnorm * sd(y)
    attr(ret, "dfs") = attr(ans$b, "df")
    attr(ret, "lambdas") = lambdas
    attr(ret, "alphas") = alphas
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
