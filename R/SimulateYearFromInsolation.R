#' @title Simulate a year from a parametric relationship to insolation
#' @param kyear point in time (kyr BP)
#' @param transfer transfer function
#' @param latitude latitude of the insolation calculation
#' @param bPolynomial (TRUE) use polynomial, FALSE = linear
#' @param anomalyRef NULL for absolute mode or vector of 365 anomaly correction values
#' @return 365 predicted values for each day of the year
#' @author Thomas Laepple
#' @examples
#' library(orbitalforcing)
#' clim <- orbitalforcing:::SelSpace3D(sat.ncep.clim, lat1 = 40, lon1 = 200)
#' transfer <- EstimateTransferfunctionInsolation(clim, 40, bPlot = FALSE, b3plot = FALSE)
#' plot(1:365, SimulateYearFromInsolation(0, transfer, 40), type = "l",
#'      main = "predicted at 0 and 10kyr BP", xlab = "day of the year")
#' lines(1:365, SimulateYearFromInsolation(10, transfer, 40), col = "red")
#' @export
SimulateYearFromInsolation <- function(kyear = 0, transfer, latitude, bPolynomial = FALSE,
                                       anomalyRef = NULL) {
  insol <- TLag(DailyInsolation(kyear, latitude, 1:365)$Fsw, transfer$lag)

  if (bPolynomial) {

    ins.1 <- insol
    ins.2 <- insol^2
    ins.3 <- insol^3
    if (length(transfer$coeff) == 4)
      result <- transfer$coeff[1] + transfer$coeff[2] * ins.1 + transfer$coeff[3] *
      ins.2 + transfer$coeff[4] * ins.3
    if (length(transfer$coeff) == 3)
      result <- transfer$coeff[1] + transfer$coeff[2] * ins.1 + transfer$coeff[3] *
      ins.2
  } else {
    result <- transfer$coeff.lin[1] + transfer$coeff.lin[2] * insol

  }
  if (!is.null(anomalyRef))
    result = result + anomalyRef
  return(result)
}


##' @title timelag with rolling over at 365 days
##' @param data vector of values at 365 days
##' @param ilag lag in days
##' @return lagged values
##' @author Thomas Laepple
TLag <- function(data, ilag) {

  temp <- rep(data, 3)

  return(temp[366:730 - ilag])

}




