\name{varCompTest.control}
\alias{varCompTest.control}
\alias{varCompTests}
\alias{varCompScoreTests}
\alias{LinScoreWeightingMethods}
\title{
	Control of variance component testing
}
\description{
\code{varCompTest.control} provide a list of arguments that control the details of the testing method for testing the nullity of variance components in \code{\link{varComp.test}}. 
}
\usage{
varCompTest.control(test = "LinScore", LinScore.wt = "InvSTD", LinScore.acc = 1e-08, 
                    LinScore.lim = 1e+06L, LinScore.method = c("AS155", "SSAS155"), 
		    VM03.method = c("SSChiBarSq", "pboot"), VM03.nsim = 10000L, 
		    SS95.method = c("pboot", "pboot"), SS95.nsim = 10000L, 
		    RLRT.method = c("exact", "pboot"), RLRT.nsim = 10000L, 
		    information = 'EI'
)
varCompTests
varCompScoreTests
LinScoreWeightingMethods
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{
	A character vector, specifying the tests to be performed. Available choices are given in \code{varCompTests}. Currently, the following tests are supported: \code{'LinScore'} (linear score tests of Qu et al. [2013+]), \code{'VM03'} (projected quadratic score test of Verbeke and Molenberghs [2003]), \code{'SS95'} (projected quadratic score test of Silvapulle and Silvapulle [1995] and Hall and Praestgaard [2001]), \code{'HP01'} (equivalent to \code{'SS95'}), \code{'RLRT'} (restricted likelihood ratio test of Crainiceanu and Ruppert [2003] or its pseudo-likelihood heuristic variant of Greven et al [2008]).  
}
  \item{LinScore.wt}{
	A character vector, specifying the weighting method for the linear score test. Available choices are given in \code{LinScoreWeightingMethods}. Currently, the following are supported: \code{'EqWt'} (equal weight for every variance component), \code{'InvSTD'} (inverse standard deviation weighting), \code{'InvSqrtV'} (inverse square root of covariance matrix) and \code{'MinVar'} (weights that minimize the variance of the resulting test statistic). See Qu et al (2013+) for detailed explanations.
}
  \item{LinScore.acc}{
	A single small positive number, specifying the accuracy when computing p-values from the linear score test. This is the same as the \code{acc} argument of \code{\link[CompQuadForm:davies]{CompQuadForm::davies}}.
}
  \item{LinScore.lim}{
	A large positive integer, specifying the maximum number of integration terms when computing the \eqn{p}-value from the linear score test. This is the same as the \code{lim} argument of \code{\link[CompQuadForm:davies]{CompQuadForm::davies}}. 
}
  \item{LinScore.method}{
	A character vector with two elements. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Available choices are \code{'AS155'} (Applied Statistics algorithm 155 of Davies (1980), i.e., \code{\link[CompQuadForm:davies]{CompQuadForm::davies}}), \code{'exact'} (synonym of \code{'AS155'}), \code{'Davies'} (another synonym of \code{'AS155'}), \code{'SSAS155'} (shifted and scaled version of \code{'AS155'}, applicable to the 2nd component only), \code{'Satterthwaite'} (scaled Chi-square approximation, applicable to the 2nd component only), and \code{'Normal'} (normal approximation, applicable to the 2nd component only). 
}
  \item{VM03.method}{
	A character vector with two components for the VM03 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Currently, the 2nd component is discarded. Available choices are \code{'ChiBarSq'} (asymptotic null distribution of Chi-bar-squares; c.f., \code{\link{pchibarsq}}) and \code{'pboot'} (parametric bootstrap, i.e., the Monte Carlo approximation to the finite-sample null distribution; c.f., Qu et al [2013+]). 
}
  \item{VM03.nsim}{
	A large positive integer, specifying the number of simulations used to approximate the null, when \code{VM03.method[1]} is set to \code{'pboot'}).
}
  \item{SS95.method}{
  A vector of two components for the SS95 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Currently available chioces are \code{"ChiBarSq"} (Chi-bar-square asymptotic null; c.f., \code{\link{pchibarsq}}) and \code{"pboot"} (parametric bootstrap, i.e., the Monte Carlo approximation to the finite-sample null distribution; c.f., Qu et al [2013+])
}
  \item{SS95.nsim}{
	A large positive integer, specifying the number of simulations used to approximate the null, when \code{SS95.method[1]} is set to \code{'pboot'}).
}
  \item{RLRT.method}{
	A vector of two components for the RLRT test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Currently available choices are "exact" (for either the exact null of Crainiceanu and Ruppert [2003] or  the pseudo-likelihood heuristic method of Greven et al [2008]) and "pboot" (for the 1st component, this is an alias of "exact"; for the 2nd component, this is not yet implemented). See also \code{\link[RLRsim:RLRTSim]{RLRsim::RLRTSim}}.
}
  \item{RLRT.nsim}{
  A large positive integer, specifying the number of simulations used to approximate the null, when \code{RLRT.method[1]} is set to \code{'exact'}).
}
%  \item{Wald.method}{
%%     ~~Describe \code{Wald.method} here~~
%}
%  \item{RWD88.method}{
%%     ~~Describe \code{RWD88.method} here~~
%}
	\item{information}{
		A character scalar, specifying the information matrix to be used in score-based tests. Available choices are given in \code{informationTypes}. By default, the expected information (\code{"EI"}) is used. C.f., Verbeke and Molenberghs (2007).
}
	
}
\value{
A list of class \code{varCompTest.control} with length equal to the length \code{test} argument. This can be passed as the \code{control} or \code{fit.control} arguments to \code{\link{varComp.test}}
}
\references{
	Qu L, Guennel T, Marshall SL. (2013+) Linear Score Tests for Variance Components in Linear Mixed Models and Applications to Genetic Association Studies. \emph{Biometrics}, under revision. 
	
	Daniel B. Hall and Jens T. Praestgaard (2001) Order-restricted tests for homogeneity in generalised linear and nonlinear mixed models. \emph{Biometrika} 88, 739--751. 
	
	Mervyn J. Silvapulle and Paramsothy Silvapulle (1995) A score test against one-sided alternatives. \emph{Journal of the American Statistical Association} 90, 342--349.
	
	Geert Verbeke and Geert Molenberghs (2003) The use of score tests for inference on variance components. \emph{Biometrics} 59, 254--262. 
	
	Geert Verbeke and Geert Molenberghs (2007) What Can Go Wrong With the Score Test?. \emph{The American Statistician} 61, 289--290. 

	Ciprian M. Crainiceanu and David Ruppert (2003) Likelihood ratio tests in linear mixed models with one variance component. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} 66, 165--185.
	
	Greven S, Crainiceanu CM Kuchenhoff H and Peters A (2008) Restricted Likelihood Ratio Testing for Zero Variance Components in Linear Mixed Models. \emph{Journal of Computational and Graphical Statistics} 17, 870--891.
	
	Robert B. Davies (1980) The Distribution of a Linear Combination of \eqn{\chi^2} Random variables. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)} 29, 323--333.
}
\author{
Long Qu
}
\seealso{
	\code{\link{varComp.test}}, \code{\link[CompQuadForm:davies]{CompQuadForm::davies}}, \code{\link[RLRsim:RLRTSim]{RLRsim::RLRTSim}}
}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
