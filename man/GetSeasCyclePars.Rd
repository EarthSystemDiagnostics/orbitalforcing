% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-seasonal-cycle-parameters.R
\name{GetSeasCyclePars}
\alias{GetSeasCyclePars}
\title{Get seasonal temperature cycle parameters using GetAmpSeasCycle}
\usage{
GetSeasCyclePars(lon, lat, age.range.kyear = c(0, 25))
}
\arguments{
\item{lon}{longitude}

\item{lat}{latitude}

\item{age.range.kyear}{age in ka of the start and end of the precessionary cycle.}
}
\value{
list of seasonal cycle parameters
}
\description{
Get seasonal temperature cycle parameters using GetAmpSeasCycle
}
\examples{
GetSeasCyclePars(lon = 120, lat = 25, age.range.kyear = c(21, 42))
}
\author{
Andrew Dolman <andrew.dolman@awi.de>
}
