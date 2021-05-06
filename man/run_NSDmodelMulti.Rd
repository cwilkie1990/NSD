\name{run.NSDmodelMulti}
\alias{run.NSDmodelMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric statistical downscaling modelling
}
\description{
This function carries out nonparametric statistical downscaling for the given inputs (with one "response" dataset and two "covariates"), resulting in simulated values from the posterior distributions of each parameter and of the resulting predictions. Typically, the "response" variable is a point-location, point-time dataset that is assumed to be accurate within measurement error, while the "covariates" are datasets of remotely-sensed grid-scale, (possibly) temporally averaged data with better spatial and/or temporal coverage than the response dataset. Examples of two covariates include satellite and drone data of the same variable. Note that this model had not yet been included in any publications.
}
\usage{
run.NSDmodelMulti(nIter = 100, nBurnIn = 10, nChains = 2, nThin = 10, yData, 
                  xData, zData, xPred, zPred, coordsData, coordsPred, 
                  times.yData = NULL, times.xData = NULL, times.zData = NULL, 
                  times.xPred = NULL, times.zPred = NULL, times.yPred = NULL, 
                  basis.type = "bspline", basis.dim = 5, period = 1, 
                  phiAlpha = 0.1, phiBeta = 0.1, phiGamma = 0.1, aAlpha = 2, 
                  bAlpha = 1, aBeta = 2, bBeta = 1, aGamma = 2, bGamma = 1, 
                  aY = 2, bY = 1, aC = 2, bC = 1, aX = 2, bX = 1, aZ = 2, 
                  bZ = 1, muD = NULL, SigmaD = NULL, muE = NULL, SigmaE = NULL, 
                  sigmaAlphaPrecInit = 1, sigmaBetaPrecInit = 1,
                  sigmaGammaPrecInit = 1, sigmaYPrecInit = 1, 
                  sigmaCPrecInit = 1, alphaInit = NULL, betaInit = NULL, 
                  gammaInit = NULL, cInit = NULL, sigmaXPrecInit = 1, 
                  sigmaZPrecInit = 1, dInit = NULL, eInit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nIter}{
Integer: the number of iterations for which to run the MCMC computations, with a default value of 100.
}
  \item{nBurnIn}{
Integer: the number of iterations to discard at the beginning of the MCMC computations, with a default value of 10.
}
  \item{nChains}{
Integer: the number of chains for which to run the MCMC computations, with a default value of 2.
}
  \item{nThin}{
Integer: the number of iterations to thin by (i.e. every (\code{nThin})th value is saved), with a default value of 10.
}
  \item{yData}{
Either: (1) Matrix of response data, of dimension \eqn{q} times by \eqn{n} data locations, with missing data set to NA; or (2) List of response data, of length equal to the number of response data locations, \eqn{n}. Each component of the list is a numeric vector of the response data at that location, of length \eqn{q_i}.
}
  \item{xData}{
Either: (1) Matrix of data for covariate 1 for the grid cells containing the response data locations, of dimension \eqn{p} times by \eqn{n} locations, with missing data set to NA; or (2) List of data for covariate 1, of length equal to the number of response data locations \eqn{n}. Each component of the list is a numeric vector of the data for covariate 1, for the grid cell containing the corresponding response data location, of length \eqn{p_i}.
}
  \item{zData}{
Either: (1) Matrix of data for covariate 2 for the grid cells containing the response data locations, of dimension \eqn{r} times by \eqn{n} locations, with missing data set to NA; or (2) List of data for covariate 2, of length equal to the number of response data locations \eqn{n}. Each component of the list is a numeric vector of the data for covariate 2, for the grid cell containing the corresponding response data location, of length \eqn{r_i}.
}
  \item{xPred}{
Either: (1) Matrix of data for covariate 1, for the grid cells containing the prediction locations, of dimension \eqn{\tilde{p}} times by \eqn{\tilde{n}} locations, with missing data set to NA; or (2) List of data for covariate 1, of length equal to the number of prediction locations \eqn{\tilde{n}}. Each component of the list is a numeric vector of the data for covariate 1, for the grid cell containing the corresponding prediction location, of length \eqn{\tilde{p}_i}.
}
  \item{zPred}{
Either: (1) Matrix of data for covariate 2, for the grid cells containing the prediction locations, of dimension \eqn{\tilde{r}} times by \eqn{\tilde{n}} locations, with missing data set to NA; or (2) List of data for covariate 2, of length equal to the number of prediction locations \eqn{\tilde{n}}. Each component of the list is a numeric vector of the data for covariate 2, for the grid cell containing the corresponding prediction location, of length \eqn{\tilde{r}_i}.
}
  \item{coordsData}{
Matrix of coordinates of the response data locations, of dimension \eqn{n} by 2
}
  \item{coordsPred}{
Matrix of coordinates of the prediction locations, of dimension \eqn{\tilde{n}} by 2
}
  \item{times.yData}{
Either (1) numeric vector of times of response data observation, of length \eqn{q} with each timepoint corresponding to a row of \eqn{q} by \eqn{n} matrix \code{yData}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{q_i} corresponding to the times of response data collection for location \eqn{i}. Note that \code{times.yData} should be a numeric vector if \code{yData} is a matrix and a list if \code{yData} is a list.
}
  \item{times.xData}{
Either (1) numeric vector of times of data for covariate 1, of length \eqn{p} with each timepoint corresponding to a row of \eqn{p} by \eqn{n} matrix \code{xData}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{p_i} corresponding to the times of data for covariate 1, for location \eqn{i}. Note that \code{times.xData} should be a numeric vector if \code{xData} is a matrix and a list if \code{xData} is a list.
}
  \item{times.zData}{
Either (1) numeric vector of times of data for covariate 2, of length \eqn{r} with each timepoint corresponding to a row of \eqn{r} by \eqn{n} matrix \code{xData}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{r_i} corresponding to the times of data for covariate 2 for location \eqn{i}. Note that \code{times.zData} should be a numeric vector if \code{zData} is a matrix and a list if \code{zData} is a list.
}
  \item{times.xPred}{
Either (1) numeric vector of times of data for covariate 1, of length \eqn{\tilde{p}} with each timepoint corresponding to a row of \eqn{\tilde{p}} by \eqn{n} matrix \code{xPred}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{\tilde{p}_i} corresponding to the times of data for covariate 1, for location \eqn{i}. Note that \code{times.xPred} should be a numeric vector if \code{xPred} is a matrix and a list if \code{xPred} is a list.
}
  \item{times.zPred}{
Either (1) numeric vector of times of data for covariate 2, of length \eqn{\tilde{r}} with each timepoint corresponding to a row of \eqn{\tilde{r}} by \eqn{n} matrix \code{zPred}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{\tilde{r}_i} corresponding to the times of data for covariate 2, for location \eqn{i}. Note that \code{times.zPred} should be a numeric vector if \code{zPred} is a matrix and a list if \code{zPred} is a list.
}
  \item{times.yPred}{
Numeric vector of times at which to predict, of length \eqn{\tilde{q}}.
}
  \item{basis.type}{
Either \code{"bspline"} or \code{"fourier"}, indicating whether to use a B-spline or Fourier basis in model fitting. See \code{\link[fda]{create.bspline.basis}} or \code{\link[fda]{create.fourier.basis}} for more details. Defaults to \code{"bspline"}.
}
  \item{basis.dim}{
Dimension of B-spline or Fourier basis used in model fitting. See \code{\link[fda]{create.bspline.basis}} or \code{\link[fda]{create.fourier.basis}} for more details. Defaults to 5.
}
  \item{period}{
Period of Fourier basis. Ignored if \code{basis.type = "bspline"}. See \code{\link[fda]{create.fourier.basis}} for more details. Defaults to 1. Note that \code{run.NSDmodelMulti} converts Dates to numeric vectors, so that a period of one year should be entered as 1 rather than as 365 or 366.
}
  \item{phiAlpha}{
Numeric: the value of the spatial decay parameter \eqn{\phi_{\alpha}}, with default value 0.1.
}
  \item{phiBeta}{
Numeric: the value of the spatial decay parameter \eqn{\phi_{\beta}}, with default value 0.1.
}
  \item{phiGamma}{
Numeric: the value of the spatial decay parameter \eqn{\phi_{\gamma}}, with default value 0.1.
}
  \item{aAlpha}{
Numeric: the value of the shape parameter \eqn{a_\alpha}, with default value 2.
}
  \item{bAlpha}{
Numeric: the value of the rate parameter \eqn{b_\alpha}, with default value 1.
}
  \item{aBeta}{
Numeric: the value of the shape parameter \eqn{a_\beta}, with default value 2.
}
  \item{bBeta}{
Numeric: the value of the rate parameter \eqn{b_\beta}, with default value 1.
}
  \item{aGamma}{
Numeric: the value of the shape parameter \eqn{a_\gamma}, with default value 2.
}
  \item{bGamma}{
Numeric: the value of the rate parameter \eqn{b_\gamma}, with default value 1.
}
  \item{aY}{
Numeric: the value of the shape parameter \eqn{a_Y}, with default value 2.
}
  \item{bY}{
Numeric: the value of the rate parameter \eqn{b_Y}, with default value 1.
}
  \item{aC}{
Numeric: the value of the shape parameter \eqn{a_c}, with default value 2.
}
  \item{bC}{
Numeric: the value of the rate parameter \eqn{b_c}, with default value 1.
}
  \item{aX}{
Numeric: the value of the shape parameter \eqn{a_x}, with default value 2.
}
  \item{bX}{
Numeric: the value of the rate parameter \eqn{b_x}, with default value 1.
}
  \item{aZ}{
Numeric: the value of the shape parameter \eqn{a_z}, with default value 2.
}
  \item{bZ}{
Numeric: the value of the rate parameter \eqn{b_z}, with default value 1.
}
  \item{muD}{
Numeric vector: the mean vector \eqn{\boldsymbol{\mu}_d}, of length \eqn{m}. If not provided, \code{muD} defaults to \code{rep(0, basis.dim)}.
}
  \item{SigmaD}{
Matrix: the covariance matrix \eqn{\boldsymbol{\Sigma}_d}, of dimension \eqn{m} by \eqn{m}. If not provided, \code{SigmaD} defaults to \code{100 * diag(basis.dim)}.
}
  \item{muE}{
Numeric vector: the mean vector \eqn{\boldsymbol{\mu}_e}, of length \eqn{m}. If not provided, \code{muE} defaults to \code{rep(0, basis.dim)}.
}
  \item{SigmaE}{
Matrix: the covariance matrix \eqn{\boldsymbol{\Sigma}_e}, of dimension \eqn{m} by \eqn{m}. If not provided, \code{SigmaE} defaults to \code{100 * diag(basis.dim)}.
}
  \item{sigmaAlphaPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{\alpha}^2)^{-1}}, with default value 1.
}
  \item{sigmaBetaPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{\beta}^2)^{-1}}, with default value 1.
}
  \item{sigmaGammaPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{\gamma}^2)^{-1}}, with default value 1.
}
  \item{sigmaYPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{y}^2)^{-1}}, with default value 1.
}
  \item{sigmaCPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{c}^2)^{-1}}, with default value 1.
}
  \item{alphaInit}{
Matrix: the initial value of \eqn{\boldsymbol{\alpha}}. If not provided, \code{alphaInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{betaInit}{
Matrix: the initial value of \eqn{\boldsymbol{\beta}}. If not provided, \code{betaInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{gammaInit}{
Matrix: the initial value of \eqn{\boldsymbol{\gamma}}. If not provided, \code{gammaInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{cInit}{
Matrix: the initial value of \eqn{\textbf{c}}. If not provided, \code{cInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{sigmaXPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{x}^2)^{-1}}, with default value 1.
}
  \item{sigmaZPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{z}^2)^{-1}}, with default value 1.
}
  \item{dInit}{
Matrix: the initial value of \eqn{\textbf{d}}. If not provided, \code{dInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{eInit}{
Matrix: the initial value of \eqn{\textbf{e}}. If not provided, \code{eInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an \code{mcmc.list} object of length \code{nChains}. Each component is an \code{mcmc} object containing the posterior simulations for each parameter in the model (including predictions).
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
## Predictions over time for location P1 (see Wilkie et al. (2019)):
if(require("fda")){
  data(ISdata)
  data(RSdata)
  data(coords.IS)
  data(coords.RS)
  data(coords.outline)
  data(which.closest)
  data(which.closest.RS)
  data(months.balaton)

  # Let's predict at a sequence of 250 times at prediction location P1:
  et.pred <- seq(min(months.balaton), max(months.balaton), length.out = 250)
  p.RS.1 <- 425
  pred.P1 <- run.NSDmodelMulti(nIter = 1000, nBurnIn = 10, nChains = 2, nThin = 2,
                   yData = ISdata, xData = RSdata[, which.closest], 
                   zData = RSdata[, which.closest], xPred = RSdata[, p.RS.1], 
                   zPred = RSdata[, p.RS.1], coordsData = coords.IS, 
                   coordsPred = coords.RS[p.RS.1, ], times.yData = months.balaton, 
                   times.xData = months.balaton, times.xPred = months.balaton, 
                   times.yPred = et.pred, times.zData = months.balaton, 
                   times.zPred = months.balaton, basis.type = "fourier",
                   basis.dim = 9, period = 1)
  # Warnings appear about providing times in Date format. This is OK, but note 
  # that we set period = 1 instead of period = 365, since the function converts 
  # the dates to decimal years, e.g. "2002-06-15" is converted to approximately 
  # 2002.452, so that one unit of time is one year instead of one day.

  summary.P1 <- summary_NSD(pred.P1)

  # Get the predictions and 95\% credible interval bounds from the summary:

  # Plot the estimates for P1 and the 95\% credible intervals:
  plot(summary.P1$pred.mat[, 1] ~ et.pred, type = "l", 
       ylim = c(min(summary.P1$lwrbnd.mat[, 1]), max(summary.P1$uprbnd.mat[, 1])), 
       xlab = "Year", ylab = "Log(chlorophyll-a)")
  lines(summary.P1$lwrbnd.mat[, 1] ~ et.pred, lty = 2)
  lines(summary.P1$uprbnd.mat[, 1] ~ et.pred, lty = 2)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }% use one of  RShowDoc("KEYWORDS")
\keyword{ IO }% __ONLY ONE__ keyword per line
