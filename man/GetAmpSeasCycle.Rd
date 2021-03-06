% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-seasonal-cycle-parameters.R
\name{GetAmpSeasCycle}
\alias{GetAmpSeasCycle}
\title{Estimate the amplitude of the seasonal temperature cycle from insolation
and a transfer function to modern surface temperatures.}
\usage{
GetAmpSeasCycle(lon, lat, kyear, transfer = NULL, plot = FALSE)
}
\arguments{
\item{lon}{longitude}

\item{lat}{latitude}

\item{kyear}{timepoint in ka}

\item{plot}{logical, plot the modern and inferred seasonal cycles}
}
\value{
numeric, amplitude of seasonal cycle
}
\description{
Estimate the amplitude of the seasonal temperature cycle from insolation
and a transfer function to modern surface temperatures.
}
\examples{
GetAmpSeasCycle(120, c(25), kyear = 5, plot = TRUE)
}
\author{
Andrew Dolman <andrew.dolman@awi.de>
}
