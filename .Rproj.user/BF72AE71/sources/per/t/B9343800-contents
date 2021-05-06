\name{NSD_data}
\alias{NSD_data}
\alias{ISdata}
\alias{RSdata}
\alias{coords.IS}
\alias{coords.RS}
\alias{coords.outline}
\alias{which.closest}
\alias{which.closest.RS}
\alias{months.balaton}
\docType{data}
\title{
Data files included in the NSD package.
}
\description{
The data files included in the NSD package allow the user to run an example analysis using nonparametric statistical downscaling.
}
\usage{
data("ISdata")
data("RSdata")
data("coords.IS")
data("coords.RS")
data("coords.outline")
data("which.closest")
data("which.closest.RS")
data("months.balaton")
}
\format{
\describe{
\item{\code{ISdata}}{\code{num [1:115, 1:9] NA NA NA NA NA NA NA NA NA NA ...}}
\item{\code{RSdata}}{\code{num [1:115, 1:7616] 1.04 3.03 5.03 5.31 4.71 ...}}
\item{\code{coords.IS}}{\code{num [1:9, 1:2] 17.3 17.4 17.7 17.9 18 ...}}
\item{\code{coords.RS}}{\code{num [1:7616, 1:2] 17.3 17.3 17.3 17.3 17.3 ...}}
\item{\code{coords.outline}}{\code{num [1:77, 1:2] 17.4 17.4 17.4 17.5 17.5 ...}}
\item{\code{which.closest}}{\code{num [1:9] 146 1005 3432 4833 6122 ...}}
\item{\code{which.closest.RS}}{\code{num [1:997] 793 1017 1095 1310 1410 ...}}
\item{\code{months.balaton}}{\code{Date[1:115], format: "2002-06-15" "2002-07-15" "2002-08-14" "2002-09-14" "2002-10-14" "2002-11-14" "2002-12-14" "2003-01-14" "2003-02-13" ...}}
}

}
\details{
\describe{
\item{\code{ISdata}}{Matrix of the in-situ log(chlorophyll-a) data for 115 months and 9 point locations in Lake Balaton, Hungary.}
\item{\code{RSdata}}{Matrix of the remote sensing log(chlorophyll-a) data for 115 months and 7616 grid cells in Lake Balaton, Hungary.}
\item{\code{coords.IS}}{Matrix of coordinates (longitude/latitude) corresponding to the 9 in-situ data point locations.}
\item{\code{coords.RS}}{Matrix of coordinates (longitude/latitude) corresponding to the 7616 remote sensing grid cell centres.}
\item{\code{coords.outline}}{Matrix of coordinates (longitude/latitude) corresponding to 77 point locations defining the lake boundaries.}
\item{\code{which.closest}}{Vector of which remote sensing grid cell is closest in geographical space to each in-situ point location.}
\item{\code{which.closest.RS}}{Vector of which remote sensing grid cell is closest in geographical space to 997 points defined by a Delaunay triangulation of the lake, constrained by the points on the lake edges (\code{coords.outline}).}
\item{\code{months.balaton}}{Date vector of 115 mid-month dates (strictly speaking equally spaced dates between \code{2002-06-15} and \code{2011-07-15}) corresponding to the data in ISdata and RSdata.}
}
}
\source{
All data in this package were made available by the GloboLakes project (\url{http://www.globolakes.ac.uk/}). The in-situ data were obtained from the Balaton Limnological Institute and the Central Transdanubian (Regional) Inspectorate for Environmental Protection, Nature Conservation and Water Management. The remote sensing data were obtained from the Diversity II project (\url{http://www.diversity2.info/products/inlandwaters/}) as part of GloboLakes. They were converted from Earth-surface reflectance observations of the MERIS instrument on board the European Space Agency's ENVISAT satellite, via algorithms.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ISdata)
data(RSdata)
data(coords.IS)
data(coords.RS)
data(coords.outline)
data(which.closest)
data(which.closest.RS)
data(months.balaton)

# Plot the in-situ and remote sensing data locations:

plot(coords.RS, xlab = "Longitude (degrees East)", ylab = "Latitude (degrees North)",
  pch = 16, cex = 0.25)
points(coords.IS, pch = 16, cex = 0.5, col = 2)
points(coords.RS[which.closest, ], pch = 2, cex = 1, col = 4)
legend("topleft", legend = c("Remote sensing grid cell centres",
  "In-situ data point locations", "Closest grid cell to each point location"),
  pch = c(16, 16, 2), pt.cex = c(0.25, 0.5, 1), col = c(1, 2, 4))

# Plot the prediction locations:

plot(coords.RS, pch = 16, col = "grey", xlab = "Longitude (degrees East)",
  ylab = "Latitude (degrees North)")
points(coords.outline, pch = 16, col = 2, cex = 0.75)
points(coords.RS[which.closest.RS, ], pch = 16, col = 1, cex = 0.5)
legend("topleft", legend = c("Remote sensing grid cells (all 7617)",
  "Constraining boundary points (chosen manually)", 
  "Remote sensing grid cells (997 selected by Delaunay triangulation)"),
  pch = 16, pt.cex = c(1, 0.75, 0.5),
  col = c("grey", "red", "black"))

# Plot the spatial patterns in log(chlorophyll-a) for one month:

if(require("sp")){
  RSdata.sp <- cbind.data.frame(coords.RS, t(RSdata))
  colnames(RSdata.sp) <- c("lon", "lat", paste0("x.", 1:115))
  sp::coordinates(RSdata.sp) <- c("lon", "lat")
  sp::spplot(RSdata.sp, "x.106", col.regions = rev(heat.colors(100)), colorkey = TRUE,
    scales = list(draw = TRUE), xlab = "Longitude (degrees East)",
    ylab = "Latitude (degrees North)")
} # Month 106 is chosen as an example.

# Plot the in-situ data for all months against the remote sensing data for their
# corresponding grid cells:

plot(ISdata, RSdata[, which.closest], xlab = "In-situ log(chlorophyll-a)", 
  ylab = "Remote sensing log(chlorophyll-a)"); abline(0,1)

}
\keyword{datasets}
