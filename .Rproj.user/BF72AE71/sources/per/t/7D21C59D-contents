\name{run.NSDmodel}
\alias{run.NSDmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric statistical downscaling modelling
}
\description{
This function carries out nonparametric statistical downscaling for the given inputs, resulting in simulated values from the posterior distributions of each parameter and of the resulting predictions. Typically, the "response" variable is a point-location, point-time dataset that is assumed to be accurate within measurement error, while the "covariate" is a dataset of remotely-sensed grid-scale, (possibly) temporally averaged dataset with better spatial and/or temporal coverage than the response dataset.
}
\usage{
run.NSDmodel(nIter = 100, nBurnIn = 10, nChains = 2, nThin = 10, yData, 
             xData, xPred, coordsData, coordsPred, times.yData = NULL, 
             times.xData = NULL, times.xPred = NULL, times.yPred = NULL, 
             basis.type = "bspline", basis.dim = 5, period = 1, 
             phiAlpha = 0.1, phiBeta = 0.1, aAlpha = 2, bAlpha = 1, 
             aBeta = 2, bBeta = 1, aY = 2, bY = 1, aC = 2, bC = 1, aX = 2,
             bX = 1, muD = NULL, SigmaD = NULL, sigmaAlphaPrecInit = 1, 
             sigmaBetaPrecInit = 1, sigmaYPrecInit = 1, sigmaCPrecInit = 1,
             alphaInit = NULL, betaInit = NULL, cInit = NULL, 
             sigmaXPrecInit = 1, dInit = NULL)
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
Either: (1) Matrix of response data, of dimension \eqn{q} times by \eqn{n} \emph{in situ} data locations, with missing data set to NA; or (2) List of response data, of length equal to the number of response data locations, \eqn{n}. Each component of the list is a numeric vector of the response data at that location, of length \eqn{q_i}.
}
  \item{xData}{
Either: (1) Matrix of covariate data for the grid cells containing the response locations, of dimension \eqn{p} times by \eqn{n} locations, with missing data set to NA; or (2) List of covariate data, of length equal to the number of response data locations \eqn{n}. Each component of the list is a numeric vector of the covariate data for the grid cell containing the corresponding response data location, of length \eqn{p_i}.
}
  \item{xPred}{
Either: (1) Matrix of covariate data for the grid cells containing the prediction locations, of dimension \eqn{\tilde{p}} times by \eqn{\tilde{n}} locations, with missing data set to NA; or (2) List of covariate data, of length equal to the number of prediction locations \eqn{\tilde{n}}. Each component of the list is a numeric vector of the covariate data for the grid cell containing the corrsponding prediction location, of length \eqn{\tilde{p}_i}.
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
Either (1) numeric vector of times of covariate data, of length \eqn{p} with each timepoint corresponding to a row of \eqn{p} by \eqn{n} matrix \code{xData}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{p_i} corresponding to the times of covariate data for location \eqn{i}. Note that \code{times.xData} should be a numeric vector if \code{xData} is a matrix and a list if \code{xData} is a list.
}
  \item{times.xPred}{
Either (1) numeric vector of times of covariate data, of length \eqn{\tilde{p}} with each timepoint corresponding to a row of \eqn{\tilde{p}} by \eqn{n} matrix \code{xPred}; or (2) List of length \eqn{n}, with each member being a numeric vector of length \eqn{\tilde{p}_i} corresponding to the times of covariate data for location \eqn{i}. Note that \code{times.xPred} should be a numeric vector if \code{xPred} is a matrix and a list if \code{xPred} is a list.
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
Period of Fourier basis. Ignored if \code{basis.type = "bspline"}. See \code{\link[fda]{create.fourier.basis}} for more details. Defaults to 1. Note that \code{run.NSDmodel} converts Dates to numeric vectors, so that a period of one year should be entered as 1 rather than as 365 or 366.
}
  \item{phiAlpha}{
Numeric: the value of the spatial decay parameter \eqn{\phi_{\alpha}}, with default value 0.1.
}
  \item{phiBeta}{
Numeric: the value of the spatial decay parameter \eqn{\phi_{\beta}}, with default value 0.1.
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
  \item{muD}{
Numeric vector: the mean vector \eqn{\boldsymbol{\mu}_d}, of length \eqn{m}. If not provided, \code{muD} defaults to \code{rep(0, basis.dim)}.
}
  \item{SigmaD}{
Matrix: the covariance matrix \eqn{\boldsymbol{\Sigma}_d}, of dimension \eqn{m} by \eqn{m}. If not provided, \code{SigmaD} defaults to \code{100 * diag(basis.dim)}.
}
  \item{sigmaAlphaPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{\alpha}^2)^{-1}}, with default value 1.
}
  \item{sigmaBetaPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{\beta}^2)^{-1}}, with default value 1.
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
  \item{cInit}{
Matrix: the initial value of \eqn{\textbf{c}}. If not provided, \code{cInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
}
  \item{sigmaXPrecInit}{
Numeric: the initial value of \eqn{(\sigma_{x}^2)^{-1}}, with default value 1.
}
  \item{dInit}{
Matrix: the initial value of \eqn{\textbf{d}}. If not provided, \code{dInit} defaults to \code{matrix(0.1, nrow = basis.dim, ncol =} \eqn{n}\code{)}.
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
  
  # Let's predict at a sequence of 250 times at prediction location P1 (see
  # Wilkie et al. 2019):
  et.pred <- seq(min(months.balaton), max(months.balaton), length.out = 250)
  p.RS.1 <- 425
  pred.P1 <- run.NSDmodel(nIter = 1000, nBurnIn = 10, nChains = 2, nThin = 2, 
                   yData = ISdata, xData = RSdata[, which.closest], 
                   xPred = RSdata[, p.RS.1], coordsData=coords.IS, 
                   coordsPred = coords.RS[p.RS.1, ], times.yData = months.balaton, 
                   times.xData = months.balaton, times.xPred = months.balaton, 
                   times.yPred = et.pred, basis.type = "fourier", basis.dim = 9,
                   period = 1)
  # Warnings appear about providing times in Date format. This is OK, but note that 
  # we set period = 1 instead of period = 365, since the function converts the dates 
  # to decimal years, e.g. "2002-06-15" is converted to approximately 2002.452, so 
  # that one unit of time is one year instead of one day.
  
  summary.P1 <- summary_NSD(pred.P1)
  
  # Get the predictions and 95\% credible interval bounds from the summary, and plot 
  # the estimates for P1 and the 95\% credible intervals:
  plot(summary.P1$pred.mat[,1] ~ et.pred, type = "l", 
    ylim = c(min(summary.P1$lwrbnd.mat[, 1]), max(summary.P1$uprbnd.mat[, 1])),
    xlab = "Year", ylab = "Log(chlorophyll-a)")
  lines(summary.P1$lwrbnd.mat[,1] ~ et.pred, lty = 2)
  lines(summary.P1$uprbnd.mat[,1] ~ et.pred, lty = 2)
}

\dontrun{
  ## Predictions across the lake, for the middle of each month:
  ## NOTE: This is very slow to run! ##
  if(require("fda")&require("sp")){
    data(ISdata)
    data(RSdata)
    data(coords.IS)
    data(coords.RS)
    data(coords.outline)
    data(which.closest)
    data(which.closest.RS)
    data(months.balaton)
    
    # Let us use a Fourier basis of dimension 9 to fit the model and predict at 997
    # locations as chosen by a constrained Delaunay triangulation (see Wilkie et al. 
    # 2019) and all 115 months in the dataset:
    
    pred.months.balaton.dec <- run.NSDmodel(nIter = 1000, nBurnIn = 10, nChains = 2, 
                                     nThin = 2, yData = ISdata, 
                                     xData = RSdata[, which.closest], 
                                     xPred = RSdata[, which.closest.RS], 
                                     coordsData = coords.IS, 
                                     coordsPred = coords.RS[which.closest.RS, ], 
                                     times.yData = months.balaton, 
                                     times.xData = months.balaton, 
                                     times.xPred = months.balaton, 
                                     times.yPred = months.balaton,
                                     basis.type = "fourier", basis.dim = 9,
                                     period = 1)
    
    summary.months.balaton.dec <- summary_NSD(pred.months.balaton.dec)
    
    # Get the matrices of predictions and 95\% credible interval bounds from 
    # the summary:
    pred.mat.all.months.balaton.dec <- summary.months.balaton.dec$pred.mat
    lwrbnd.mat.all.months.balaton.dec <- summary.months.balaton.dec$lwrbnd.mat
    uprbnd.mat.all.months.balaton.dec <- summary.months.balaton.dec$uprbnd.mat
    
    # Create a SpatialPointsDataFrame for the predictions, lower and upper 
    # 95\% credible interval bounds (using package "sp"):
    pred.mat.all.months.balaton.dec.sp <- 
      cbind.data.frame(coords.RS[which.closest.RS, ], 
        t(pred.mat.all.months.balaton.dec))
    colnames(pred.mat.all.months.balaton.dec.sp) <- c("lon", "lat", 
      paste0("x.", 1:length(months.balaton)))
    sp::coordinates(pred.mat.all.months.balaton.dec.sp) <- c("lon","lat")
    
    lwrbnd.mat.all.months.balaton.dec.sp <- 
      cbind.data.frame(coords.RS[which.closest.RS, ], 
        t(lwrbnd.mat.all.months.balaton.dec))
    colnames(lwrbnd.mat.all.months.balaton.dec.sp) <- c("lon", "lat",
      paste0("x.", 1:length(months.balaton)))
    sp::coordinates(lwrbnd.mat.all.months.balaton.dec.sp) <- c("lon", "lat")
    
    uprbnd.mat.all.months.balaton.dec.sp <- 
      cbind.data.frame(coords.RS[which.closest.RS, ], 
        t(uprbnd.mat.all.months.balaton.dec))
    colnames(uprbnd.mat.all.months.balaton.dec.sp) <- c("lon", "lat", 
      paste0("x.", 1:length(months.balaton)))
    sp::coordinates(uprbnd.mat.all.months.balaton.dec.sp) <- c("lon", "lat")
    
    # Plot the predictions for a few months (arranged in same figure):
    sp1 <- sp::spplot(pred.mat.all.months.balaton.dec.sp, "x.1", 
                      col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                      main = "15th June 2002", scales = list(draw = TRUE), 
                      xlab = "Longitude (degrees East)", 
                      ylab = "Latitude (degrees North)")
    sp2 <- sp::spplot(pred.mat.all.months.balaton.dec.sp, "x.2", 
                      col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                      main = "15th July 2002", scales = list(draw = TRUE), 
                      xlab = "Longitude (degrees East)", 
                      ylab = "Latitude (degrees North)")
    sp3 <- sp::spplot(pred.mat.all.months.balaton.dec.sp, "x.3", 
                      col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                      main = "14th August 2002", scales = list(draw = TRUE), 
                      xlab = "Longitude (degrees East)",
                      ylab = "Latitude (degrees North)")
    
    plot(sp1, position = c(0, 0.66, 1, 0.99), more = TRUE)
    plot(sp2, position = c(0, 0.33, 1, 0.66), more = TRUE)
    plot(sp3, position = c(0, 0, 1, 0.33), more = FALSE)
    
    # Plot the predictions for 15th June 2002, along with lower and upper 
    # 95\% credible interval bounds (on same scale):
    cuts <- seq(min(lwrbnd.mat.all.months.balaton.dec.sp@data$"x.1"), 
                max(uprbnd.mat.all.months.balaton.dec.sp@data$"x.1"), length.out = 100)
    sp1a <- sp::spplot(pred.mat.all.months.balaton.dec.sp, "x.1", 
                col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                main = "Predictions", scales = list(draw = TRUE), 
                xlab = "Longitude (degrees East)", 
                ylab = "Latitude (degrees North)", cuts = cuts)
    sp4 <- sp::spplot(lwrbnd.mat.all.months.balaton.dec.sp, "x.1", 
                col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                main = "Lower 95 CI bound", scales = list(draw = TRUE), 
                xlab = "Longitude (degrees East)", 
                ylab = "Latitude (degrees North)", cuts = cuts)
    sp5 <- sp::spplot(uprbnd.mat.all.months.balaton.dec.sp, "x.1", 
                col.regions = rev(heat.colors(100)), colorkey = TRUE, 
                main = "Upper 95 CI bound", scales = list(draw = TRUE), 
                xlab = "Longitude (degrees East)", 
                ylab = "Latitude (degrees North)", cuts = cuts)  
    plot(sp4, position = c(0, 0.66, 1, 0.99), more = TRUE)
    plot(sp1a, position = c(0, 0.33, 1, 0.66), more = TRUE)
    plot(sp5, position = c(0, 0, 1, 0.33), more = FALSE)
  }}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }% use one of  RShowDoc("KEYWORDS")
\keyword{ IO }% __ONLY ONE__ keyword per line
