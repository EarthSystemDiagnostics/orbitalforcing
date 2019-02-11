## Functions to get amplitude of seasonal cycle from insolation and transfer
## function to modern climate

#' Estimate the amplitude of the seasonal temperature cycle from insolation
#' and a transfer function to modern surface temperatures.
#'
#' @param lon longitude
#' @param lat latitude
#' @param kyear timepoint in ka
#' @param plot logical, plot the modern and inferred seasonal cycles
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @return numeric, amplitude of seasonal cycle
#' @export
#'
#' @examples
#' GetAmpSeasCycle(120, c(25), kyear = 5, plot = TRUE)
GetAmpSeasCycle <- function(lon, lat, kyear, transfer = NULL, plot = FALSE){

  stopifnot(length(lat) == 1, length(lon) == 1, length(kyear) == 1)



  if (is.null(transfer)){
    #First extract the modern daily temperature from NCEP reanalysis
    climate.modern <- ecustools::SelSpace3D(ecustools::sat.ncep.clim,
                                            lat1 = lat, lon1 = lon) - 273.15

    # Estimate the modern relationship between insolation and T
    transfer <- EstimateTransferfunctionInsolation(
      climate = climate.modern, latitude = lat)
  }


  climate.at.kyear <- PaleoSpec::MonthlyFromDaily(
    ecustools::SimulateYearFromInsolation(kyear = kyear, transfer = transfer,
                                          latitude = lat, bPolynomial = TRUE)
  )

  if (plot){
    ylims <- range(c(range(climate.modern), range(climate.at.kyear)))
    pars.mfrow <- par()$mfrow
    par(mfrow = c(1, 2))
    plot(climate.modern, ylim = ylims, xlab = "DOY")
    plot(climate.at.kyear, type = "b", ylim = ylims, xlab = "Month")
    par(mfrow = pars.mfrow)
  }

  return(diff(range(climate.at.kyear)))
}



#' Get seasonal temperature cycle parameters using GetAmpSeasCycle
#'
#' @param age.range.kyear age in ka of the start and end of the precessionary cycle
#' @inheritParams GetAmpSeasCycle
#' @return list of seasonal cycle parameters
#' @export
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' GetSeasCyclePars(lon = 120, lat = 25, age.range.kyear = c(21, 42))
GetSeasCyclePars <- function(lon, lat, age.range.kyear = c(0, 25)){

  if (diff(age.range.kyear) > 25) warning("Searching over multiple precessionary cycles")

  # Estimate the modern relationship between insolation and T
  #climate.modern <- ecustools::SelSpace3D(ecustools::sat.ncep.clim,
  #                                        lat1 = lat, lon1 = lon) - 273.15

  # transfer <- EstimateTransferfunctionInsolation(
  #   climate = climate.modern, latitude = lat)

  glob.tf.sat.ncep.clim <- glob.tf.sat.ncep.clim

  WhichTF <- function(lon, lat){

    coords <- attr(glob.tf.sat.ncep.clim, "split_labels")

    lats <- unique(coords$lat)
    lons <- unique(coords$lon)

    nearest.lat <- lats[which.min(abs(lats - lat))]
    nearest.lon <- lons[which.min(abs(lons - lon))]

    ind <- which(coords$lon == nearest.lon & coords$lat == nearest.lat)

    glob.tf.sat.ncep.clim[[ind]][[1]]

  }

  transfer <- WhichTF(lon = lon, lat = lat)

  f <- function(kyear, lon, lat) {
    GetAmpSeasCycle(lon = lon, lat = lat, kyear = kyear, transfer = transfer)
  }

  max.amp <- optimize(f, interval = age.range.kyear, lat = lat, lon = lon, maximum = TRUE)
  min.amp <- optimize(f, interval = age.range.kyear, lat = lat, lon = lon, maximum = FALSE)

  names(max.amp) <- NULL
  names(min.amp) <- NULL


  pars <- c(max.kyear = max.amp[1], max.amp = max.amp[2], min.kyear = min.amp[1], min.amp = min.amp[2])

  pars$mean.amp <- mean(c(pars$max.amp, pars$min.amp))

  pars$amp.diff <- pars$max.amp - pars$min.amp

  pars$sig.sq_c <- spectraluncertainty::VarSine(pars$mean.amp)
  pars$sig.sq_a <- ((pars$max.amp - pars$mean.amp) / pars$mean.amp)^2

  return(pars)
}


