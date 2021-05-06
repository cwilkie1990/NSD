\name{NSD-package}
\alias{NSD-package}
\docType{package}
\title{
\packageTitle{NSD}
}
\description{
\packageDescription{NSD}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{NSD}
\packageIndices{NSD}

The main function in this package is \code{\link{run.NSDmodel}}, which allows the user to run the nonparametric statistical downscaling model described in Wilkie et al. (2019). Two examples are given, showing how the function can be used to predict over both space and time, fusing data of different spatiotemporal support.

The function \code{\link{run.NSDmodelMulti}} allows fitting a nonparametric statistical downscaling model with two covariates. One example is given.

The function \code{\link{summary_NSD}} provides a list of summary statistics for each parameter in the model, along with matrices of predictions and corresponding 95\% credible interval bounds.

The package contains some example datasets for Lake Balaton, Hungary. These are described in \code{\link{NSD_data}}.

Helper functions described in \code{\link{NSD_other}} allow the user to manipulate their data into the list format required by \code{\link{run.NSDmodel}} and \code{\link{run.NSDmodelMulti}}, although these model-fitting functions can do this internally, so that it is not required to have data in list format before modelling.

%The functions described in \code{\link{NSD_helper_functions}} are used internally by \code{\link{run.NSDmodel}} and should not be called directly by the user.
}
\author{
\packageAuthor{NSD}

Maintainer: \packageMaintainer{NSD}
}
\references{
Wilkie, C.J., C.A. Miller, E.M. Scott, R.A. O'Donnell, P.D. Hunter, E. Spyrakos and A.N. Tyler (2019). \emph{Nonparametric statistical downscaling for the fusion of data of different spatiotemporal support.} Environmetrics 30(3):e2549.
}
\keyword{ package }
\seealso{
\code{\link[fda]{create.fourier.basis}}, \code{\link[fda]{create.bspline.basis}} and \code{\link[fda]{eval.basis}} are useful for creating and evaluating Fourier or B-spline bases. These functions are used within \code{\link{run.NSDmodel}} and \code{\link{run.NSDmodelMulti}}.

\code{\link[lubridate]{decimal_date}} is used within \code{\link{run.NSDmodel}} and \code{\link{run.NSDmodelMulti}} to change objects of class \code{"POSIXt"}, \code{"POSIXct"}, \code{"POSIXlt"} or \code{"Date"} to numeric vectors.

\code{\link[sp]{spplot}} provides a useful plotting interface for spatial plots (e.g. of data or of predictions).
}
\examples{
## Data manipulation:

data(ISdata)
data(months.balaton)

IS.avail<-matrix.to.list(ISdata)
times.avail<-times.to.list(ISdata,months.balaton)

IS.avail
times.avail

# Plot the spatial patterns in log(chlorophyll-a) for month 106 of 115:

data(coords.RS)
data(RSdata)

if(require("sp")){
  RSdata.sp <- cbind.data.frame(coords.RS, t(RSdata))
  colnames(RSdata.sp) <- c("lon", "lat", paste0("x.", 1:115))
  sp::coordinates(RSdata.sp) <- c("lon", "lat")
  sp::spplot(RSdata.sp, "x.106", col.regions = rev(heat.colors(100)), colorkey = TRUE,
    scales = list(draw = TRUE), xlab = "Longitude (degrees East)",
    ylab = "Latitude (degrees North)")
}

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

}
