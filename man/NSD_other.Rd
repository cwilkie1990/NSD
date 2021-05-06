\name{NSD_other}
\alias{NSD_other}
\alias{matrix.to.list}
\alias{times.to.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper functions for data manipulation, which may be called by the user.
}
\description{
These helper functions may be used to change a matrix of data into list format. These are used within the \code{run.NSDmodel} and \code{run.NSDmodelMulti} functions and do not need to be called directly.
}
\usage{
matrix.to.list(in.mat)
times.to.list(in.mat,in.times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{in.mat}{
Numeric matrix with columns representing locations and rows representing times. The matrix may contain missing values.
}
  \item{in.times}{
  Vector of times of data, corresponding to the rows of \code{in.mat}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{matrix.to.list}{This function results in a list of length equal to the number of columns in \code{in.mat} (i.e. the number of data locations). Each component of the list contains the non-missing data for that location.}
\item{times.to.list}{This function results in a list of length equal to the number of columns in \code{in.mat} (i.e. the number of data locations). Each component of the list contains the times of data collection for the data in the corresponding component of \code{matrix.to.list}.}
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
data(ISdata)
data(months.balaton)

IS.avail<-matrix.to.list(ISdata)
times.avail<-times.to.list(ISdata,months.balaton)

IS.avail
times.avail
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
