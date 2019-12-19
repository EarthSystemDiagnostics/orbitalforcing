## Functions to get amplitude of seasonal cycle from insolation and transfer
## function to modern climate


# SelSpace3D functions here as temporary dependency fix --------

#' Get name of proxy object
#'
#' The function returns the \code{"name"} attribute of a \code{"pTs"} or
#' \code{"pField"} object.
#' @param data a \code{"pTs"} or \code{"pField"} object.
#' @return Character string with the \code{"name"} attribute of \code{data}.
#' @source Function copied from "basis.R" in paleolibary/src/.
#' @author Thomas Laepple
#' @examples
#' @aliases getname GetName
#' @keywords internal
GetName <- function(data) {
  return(attr(data, "name"))
}
getname <- function(...) {
  warning("getname is deprecated and replaced with GetName
            to comply with ECUS R style guide.")
  GetName(...)
}

#' Get history of a proxy object
#'
#' The function returns the \code{"history"} attribute of a \code{"pTs"} or
#' \code{"pField"} object.
#' @inheritParams GetName
#' @family history
#' @return Character string with the \code{"history"} attribute of \code{data}.
#' @source Function copied from "basis.R" in paleolibary/src/.
#' @author Thomas Laepple
#' @examples
#' @aliases GetHistory gethistory
#' @keywords internal
GetHistory <- function(data) {
  #print(match.call())
  return(attr(data, "history"))
}
gethistory <- function(...){
  warning("gethistory is deprecated and replaced with GetHistory
            to comply with ECUS R style guide.")
  GetName(...)
}

#' Amend the history of a proxy object
#'
#' This function amends the history information of a \code{"pTs"} or
#' \code{"pField"} object by adding a given character string (e.g., a comment)
#' to the \code{"history"} attribute of the object together with the date
#' information of this action.
#' @param x a \code{"pTs"} or \code{"pField"} object.
#' @param newhist character string with the history information to be added.
#' @inheritParams pTs
#' @family history
#' @return the \code{"pTs"} or \code{"pField"} object with the amended
#' history.
#' @source Function copied from "basis.R" in paleolibary/src/.
#' @author Thomas Laepple
#' @examples
#' @aliases addhistory AddHistory
#' @keywords internal
AddHistory <- function(x, newhist, date = TRUE) {
  if (date) newhist <- paste(date(), newhist, sep = ": ")
  attr(x, "history") <- c(attr(x, "history"), newhist)
  return(x)
}
addhistory <- function(...) {
  warning("addhistory is deprecated and replaced with AddHistory
            to comply with ECUS R style guide.")
  AddHistory(...)
}


#' Create a pTs object
#'
#' @param data a vector, matrix or ts object
#' @param time vector
#' @param lat vector of length 1
#' @param lon vector of length 1
#' @param name character string, name of the pTs object to be generated
#' @param history character string to append to the history attribute
#' @param date logical, whether or not to append the current date to the history attribute
#'
#' @description pTs objects are enhanced time-series \code{\link[stats]{ts}} objects. \code{pTs()} adds attributes such as longitude and latitude to a time series
#'   vector/time series vectors (having the same time basis) and assigns the
#'   class "pTs" to the resulting object
#' @source Function copied from "proxytype.R" in paleolibary/src/
#'
#' @return pTs object
#' @keywords internal
#'
#' @author Thomas Laepple
#' @examples
pTs <- function(data,
                time,
                lat = 0,
                lon = 0,
                name = "",
                history = "",
                date = TRUE)
{
  #constants
  TOL = 1 / 400 #tolerance less than 1 day

  if (length(data) <= 1)
    if (is.null(data[1]))
      data <- matrix(NA, length(time), length(name))
    else if (length(name) == 1)
      data <- rep(data[1], length(time))
    else
      data <- matrix(data[1], length(time), length(name))



    if (!is.null(ncol(data)) && (ncol(data) > 1))
      #multiple datasets
    {
      #check data
      if (nrow(data) != length(time))
        stop("nTime != N_Elements(data)")

    }
    #check data
    else if (length(data) != (length(time)))
      stop("nTime != N_Elements(data)")

    #shape data in 2D array

    if (length(time) > 1)
      #real time series
    {
      if (abs(max(diff(time)) - min(diff(time))) > TOL)
        stop("time steps are not equidistant")
      result <- ts(data, start = time[1], deltat = (time[2] - time[1]))
    }
    else
      result <- ts(data, start = time[1]) #or only one time step

    #put attributes and classes
    attr(result, 'lat') <- lat
    attr(result, 'lon') <- lon
    attr(result, 'name') <- name

    if (date)
      attr(result, 'history') <- paste(date(), history)
    else
      attr(result, 'history') <- paste(history)


    attr(result, 'oclass') <- class(result)
    attr(result, 'nclass') <- c('pTs', 'ts')
    class(result) <- attr(result, 'nclass')

    invisible(result)
}

#' @title first element of a vector
#' @param x vector
#' @return first element of X
#' @author Thomas Laepple
#' @keywords internal
FirstElement <- function(x) {
  return(x[1])
}

#' @title last element of a vector
#' @param x vector
#' @return last element of X
#' @author Thomas Laepple
#' @keywords internal
LastElement <- function(x) {
  return(x[length(x)])
}


#' SelSpace3D
#' @param data pField object
#' @param lat1 vector length 1
#' @param lon1 vector legnth 1
#' @param SBOX vector length 1, defaults to 5
#' @param tolLon vector length 1, defaults to 10
#' @param bNN logical, defaults to FALSE
#' @param timeindexNA at which timestep should the Next Neighbour algorithm
#'   search for missing values
#'
#' @description Function to interpolate a field to a given point specified by
#'   latitude and longitude. It uses the nearest neighbour if the adjancents
#'   points are missing, if not it uses bilinear interpolation
#'
#' @return \code{\link{pTs}} object
#' @aliases SelSpace3D selspace.3D
#' @author Thomas Laepple
#' @examples
#' @keywords internal
SelSpace3D <- function(data, lat1, lon1, SBOX = 5, tolLon = 10,
                       bNN = FALSE, timeindexNA = 1) {
  choice.lat <- lat1
  choice.lon <- lon1
  temp <- attributes(data)

  if (prod(dim(data)) != length(temp$lon) * length(temp$lat) *
      length(time(data)))
    stop("N(data) != N(lat)*N(lon)*N(time)")

  # make a 2D array [lon,lat] containing the orginal indices
  pointer3d.raw <- array(1:dim(data)[2], c(1, length(temp$lon),
                                           length(temp$lat)))

  # arrange to get a continous field
  nlon <- length(temp$lon)
  nTime <- length(time(data))
  d <- diff(temp$lon)

  if (max(d) > (min(d) + 0.01)) {
    edgelon <- which(d == max(d))
    pointer3d <- array(NA, c(1, length(temp$lon), length(temp$lat)))
    pointer3d[, (edgelon + 1):nlon, ] <- pointer3d.raw[,
                                                       (edgelon + 1):nlon, ]
    pointer3d[, 1:edgelon, ] <- data3d.raw[, 1:edgelon, ]
    temp$lon <- c(temp$lon[(edgelon + 1):nlon], temp$lon[1:edgelon] +
                    360)
  } else pointer3d <- pointer3d.raw

  # longitude jump by wrapping
  wrap.dLon <- FirstElement(temp$lon) - (LastElement(temp$lon) -
                                           360)

  if ((wrap.dLon > 0) & (wrap.dLon < tolLon)) {
    ### Check if we the data is global on the longitudes, than Copy
    ### the data 3 times.... to avoid breaks in the longitude
    pointer3d.3c <- array(NA, c(1, 3 * length(temp$lon),
                                length(temp$lat)))
    pointer3d.3c[, 1:nlon, ] <- pointer3d[, 1:nlon, ]
    pointer3d.3c[, (nlon + 1):(2 * nlon), ] <- pointer3d[,
                                                         1:nlon, ]
    pointer3d.3c[, (2 * nlon + 1):(3 * nlon), ] <- pointer3d[,
                                                             1:nlon, ]

    lon.3c <- c(temp$lon[1:nlon] - 360, temp$lon[1:nlon],
                temp$lon[1:nlon] + 360)
  } else {
    # Do wrap the longitudes at the longitudes could not be
    # connected
    lon.3c <- temp$lon
    pointer3d.3c <- pointer3d
  }

  if (temp$lat[2] < temp$lat[1])
    # if the latitudes are from + to -, reverse them
  {
    temp$lat <- rev(temp$lat)
    pointer3d.3c <- pointer3d.3c[, , rev(seq(temp$lat))]
  }


  # attention... midpoints are given...
  if ((lat1 > LastElement(temp$lat)) | (lat1 < FirstElement(temp$lat))) {
    warning("Latitude outside field")
    return(NULL)
  }
  if ((lon1 > LastElement(lon.3c)) | (lon1 < FirstElement(lon.3c))) {
    warning("Longitude outside field")
    return(NULL)
  }

  indexLat <- rev(which(temp$lat <= lat1))[1]:which(temp$lat >=
                                                      lat1)[1]
  indexLon <- rev(which(lon.3c <= lon1))[1]:which(lon.3c >=
                                                    lon1)[1]

  # First check if we have any missing data
  pointer.select <- pointer3d.3c[, indexLon, indexLat]


  if (is.null(dim(pointer.select))) {
    dim(pointer.select) <- c(1, length(pointer.select))
  } else {
    dim(pointer.select) <- c(1, dim(pointer.select))
  }


  # data.select<-(as.matrix(data)[,pointer.select]) if
  # (is.null(dim(data.select)) nData<-sum(!is.na(data.select))
  # else nData<-colSums(!is.na(data.select))

  # if (min(colSums(!is.na(data)[,pointer.select]==0) {
  # warning('One neighbour point only contains missing values;
  # therefore Next Neighbour is used') bNN=TRUE }

  if (bNN) {
    # Create an area around the point
    tempindexLat <- ((indexLat[1] - SBOX):(indexLat[1] +
                                             SBOX))
    tempindexLon <- ((indexLon[1] - SBOX):(indexLon[1] +
                                             SBOX))

    # remove areas outside the boundaries
    tempindexLat <- tempindexLat[tempindexLat > 0]
    tempindexLon <- tempindexLon[tempindexLon > 0]
    tempindexLat <- tempindexLat[tempindexLat <= length(temp$lat)]
    tempindexLon <- tempindexLon[tempindexLon <= length(lon.3c)]

    # Get the nearest neighbours
    res <- expand.grid(tempindexLon, tempindexLat)

    # only retain the nearest nonmissing neighbours
    IndexRegion <- diag(pointer3d.3c[1, res[, 1], res[, 2]])
    indexV <- !is.na(data[timeindexNA, IndexRegion])

    if (sum(indexV) == 0)
      return(NA)
    x <- res[indexV, 1]
    y <- res[indexV, 2]

    # Distances to the points in the area
    D2i <- (lon.3c[x] - lon1)^2 + (temp$lat[y] - lat1)^2

    # nearest neighbour is the one with the smallest distance
    neighbour <- order(D2i)[1]

    # Switch to real data
    intpoldata <- data[, pointer3d.3c[, x[neighbour], y[neighbour]]]
    choice.lat <- temp$lat[y[neighbour]]
    choice.lon <- lon.3c[x[neighbour]]
  } else {
    # here the interpolation starts

    if (length(indexLat) == 1)
      ey = NA else ey <- (lat1 - temp$lat[indexLat[1]])/(temp$lat[indexLat[2]] -
                                                           temp$lat[indexLat[1]])
      if (length(indexLon) == 1)
        ex = NA else ex <- (lon1 - lon.3c[indexLon[1]])/(lon.3c[indexLon[2]] -
                                                           lon.3c[indexLon[1]])


        if ((!is.finite(ex)) & (!is.finite(ey)))
          intpoldata <- data[, pointer.select] else # we are on the point, no interpolation
            if (!is.finite(ex)) {
              intpoldata <- data[, pointer.select[, 1]] + (data[,
                                                                pointer.select[, 2]] - data[, pointer.select[,
                                                                                                             1]]) * ey
              # only latitudonal interpolation
            } else {
              if (!is.finite(ey)) {
                intpoldata <- data[, pointer.select[, 1]] + (data[,
                                                                  pointer.select[, 2]] - data[, pointer.select[,
                                                                                                               1]]) * ex
                # only longitudonal interpolation #Switch to real data
              } else intpoldata <- (1 - ex) * (1 - ey) * data[, pointer.select[,
                                                                               1, 1]] + (1 - ex) * (ey) * data[, pointer.select[,
                                                                                                                                1, 2]] + (ex) * (1 - ey) * data[, pointer.select[,
                                                                                                                                                                                 2, 1]] + (ex * ey) * data[, pointer.select[,
                                                                                                                                                                                                                            2, 2]]
            }
  }
  # create time series
  result <- pTs(intpoldata, time(data), choice.lat, choice.lon,
                GetName(data), GetHistory(data), date = FALSE)

  hist <- paste("selspace: lat=", lat1, " lon=", lon1, sep = "")

  return(AddHistory(result, hist))
}

#' @keywords internal
selspace.3D <- function(...) {
  warning("selspace.3D is deprecated and replaced with SelSpace3D
    to comply with ECUS R style guide")
  SelSpace3D(...)
}

# Seas cycle functions -------

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
    climate.modern <- SelSpace3D(sat.ncep.clim,
                                            lat1 = lat, lon1 = lon) - 273.15

    # Estimate the modern relationship between insolation and T
    transfer <- EstimateTransferfunctionInsolation(
      climate = climate.modern, latitude = lat)
  }


  climate.at.kyear <- PaleoSpec::MonthlyFromDaily(
    SimulateYearFromInsolation(kyear = kyear, transfer = transfer,
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
#' @param age.range.kyear age in ka of the start and end of the precessionary cycle.
#' @inheritParams GetAmpSeasCycle
#' @return list of seasonal cycle parameters
#' @export
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' GetSeasCyclePars(lon = 120, lat = 25, age.range.kyear = c(21, 42))
GetSeasCyclePars <- function(lon, lat, age.range.kyear = c(0, 25)){

  if (diff(age.range.kyear) > 25) warning("Searching over multiple precessionary cycles")

  VarSine <- function(full.amp) {
    (full.amp / 2)^2 / 2
  }


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

  pars$sig.sq_c <- VarSine(pars$mean.amp)
  pars$sig.sq_a <- ((pars$max.amp - pars$mean.amp) / pars$mean.amp)^2

  return(pars)
}


