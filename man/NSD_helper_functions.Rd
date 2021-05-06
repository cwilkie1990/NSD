\name{NSD_helper_functions}
\alias{NSD_helper_functions}
\alias{choleskyEfficient}
\alias{calcDistsArma}
\alias{H0functionArma}
\alias{NSDmodel}
\alias{NSDmodelMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper functions for the NSD package.
}
\description{
These functions are called by the function \code{run.NSDmodel} and should not be called directly by the user.
}
\usage{
choleskyEfficient(Amat)
calcDistsArma(coordsArma)
H0functionArma(phi, dists)
NSDmodel(nIter, nThin, yData, xData, xPred, coordsData, coordsPred, By, Bx, ByPred, BxPred,
         phiAlpha, phiBeta, aAlpha, bAlpha, aBeta, bBeta, aY, bY, aC, bC, aX, bX, 
         muD, SigmaD, sigmaAlphaPrecInit, sigmaBetaPrecInit, sigmaYPrecInit, 
         sigmaCPrecInit, alphaInit, betaInit, cInit, sigmaXPrecInit, dInit)
NSDmodelMulti(nIter, nThin, yData, xData, zData, xPred, zPred, coordsData, coordsPred, By, Bx, Bz, ByPred, BxPred, BzPred,
         phiAlpha, phiBeta, phiGamma, aAlpha, bAlpha, aBeta, bBeta, aGamma, bGamma, aY, bY, aC, bC, aX, bX, aZ, bZ, 
         muD, SigmaD, muE, SigmaE, sigmaAlphaPrecInit, sigmaBetaPrecInit, sigmaGammaPrecInit, sigmaYPrecInit, 
         sigmaCPrecInit, alphaInit, betaInit, gammaInit, cInit, sigmaXPrecInit, sigmaZPrecInit, dInit, eInit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\details{
\describe{
\item{\code{choleskyEfficient}}{carries out an efficient Cholesky decomposition of the matrix \code{Amat}.}
\item{\code{calcDistsArma}}{computes the Euclidean distances between the coordinates in \code{coordsArma}.}
\item{\code{H0functionArma}}{computes the exponential spatial correlation: exp(-\code{phi} \eqn{\times} \code{dists}).}
\item{\code{NSDmodel}}{carries out the MCMC computations using \code{C++}.}
\item{\code{NSDmodelMulti}}{carries out the MCMC computations using \code{C++}.}
}
}
\value{
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Craig Wilkie
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
