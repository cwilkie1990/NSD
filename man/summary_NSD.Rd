\name{summary_NSD}
\alias{summary_NSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary function for NSDmodel output
}
\description{
This function produces a summary of the output of a model fitted using the \code{run.NSDmodel} function, along with \eqn{\tilde{q}} by \eqn{\tilde{n}} matrices of predictions and corresponding 95\% credible interval bounds.
}
\usage{
summary_NSD(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class 'mcmc.list'.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of length 4 containing:
\item{summary.mcmc.list}{Usual output of coda::summary() for object of class 'mcmc.list'}
\item{pred.mat}{\eqn{\tilde{q}} by \eqn{\tilde{n}} matrix of predictions for the \eqn{\tilde{q}} prediction times and the \eqn{\tilde{n}} prediction locations}
\item{lwrbnd.mat}{\eqn{\tilde{q}} by \eqn{\tilde{n}} matrix of 95\% credible interval lower bounds for the \eqn{\tilde{q}} prediction times and the \eqn{\tilde{n}} prediction locations}
\item{uprbnd.mat}{\eqn{\tilde{q}} by \eqn{\tilde{n}} matrix of 95\% credible interval upper bounds for the \eqn{\tilde{q}} prediction times and the \eqn{\tilde{n}} prediction locations}
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
                   xPred = RSdata[, p.RS.1], coordsData = coords.IS, 
                   coordsPred = coords.RS[p.RS.1, ], times.yData = months.balaton, 
                   times.xData = months.balaton, times.xPred = months.balaton, 
                   times.yPred = et.pred, basis.type = "fourier", basis.dim = 9,
                   period = 1)
  # Warnings appear about providing times in Date format. This is OK, but note 
  # that we set period = 1 instead of period = 365, since the function converts
  # the dates to decimal years, e.g. "2002-06-15" is converted to approximately
  # 2002.452, so that one unit of time is one year instead of one day.
  
  summary.P1 <- summary_NSD(pred.P1)
  
  # Plot the estimates for P1 and the 95\% credible intervals:
  plot(summary.P1$pred.mat[, 1] ~ et.pred, type = "l", ylim = c(min(summary.P1$lwrbnd.mat[, 1]), 
       max(summary.P1$uprbnd.mat[, 1])), xlab = "Year", ylab = "Log(chlorophyll-a)")
  lines(summary.P1$lwrbnd.mat[, 1] ~ et.pred, lty = 2)
  lines(summary.P1$uprbnd.mat[, 1] ~ et.pred, lty = 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
