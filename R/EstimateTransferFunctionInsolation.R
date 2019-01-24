#' @title Estimate linear or polynomial transfer function between insolation and
#'   a climate parameter
#' @param latitude Latitude (deg North)
#' @param polynom.order order of the polynom, 2 or 3
#' @param clab name and unit = ylab of climate parameter for plotting
#' @param bPlot plot diagnostics (if TRUE)
#' @param b3plot plot more diagnostics (if TRUE)
#' @param main Title of the plot
#' @reference Laepple, T. and Lohmann, G.: Seasonal cycle as template for
#'   climate variability on astronomical timescales, Paleoceanography, 24(4),
#'   PA4201, doi:10.1029/2008PA001674, 2009.
#' @return list(coeff, coeff.lin, rsq = c(rsq1, rsq2),lag)
#' coeff coefficient of polynomial,
#' coeff.lin = coefficient of linear fit,
#' rsq of linear and of polynomial fit,
#' timelag in days
#' @author Thomas Laepple
#' @examples
#' clim <- SelSpace3D(sat.ncep.clim, lat1 = 40, lon1 = 200)
#'
#' result <- EstimateTransferfunctionInsolation(clim, 70, bPlot=TRUE, b3plot=TRUE)
#' @export
EstimateTransferfunctionInsolation <- function(climate, latitude,
  polynom.order = 3, clab = "", bPlot = FALSE, b3plot = FALSE,
  main = "") {

  climate.matrix <- matrix(NA, 365, 365)

  # Insolation -------
  insol <- DailyInsolation(0, latitude, 1:365)$Fsw

  ins.1 <- insol
  ins.2 <- insol^2
  ins.3 <- insol^3


  basis <- rep(climate, 2)

  # save lagged versions of the climate vector

  for (iopt in 1:365) climate.matrix[, iopt] <- basis[iopt:(iopt + 364)]
  if (polynom.order == 3) {
    t <- lm(climate.matrix ~ ins.1 + ins.2 + ins.3)

    # ul, il und index.neg wurde verwendet um eine positive
    # Steigung zu erreichen, nur in den Tropen wichtig

    ul <- max(insol) * t$coefficients[2, ] + max(insol)^2 *
      t$coefficients[3, ] + max(insol)^3 * t$coefficients[4, ]
    ll <- min(insol) * t$coefficients[2, ] + min(insol)^2 *
      t$coefficients[3, ] + min(insol)^3 * t$coefficients[4, ]

  } else {
    t <- lm(climate.matrix ~ ins.1 + ins.2)
    ul <- max(insol) * t$coefficients[2, ] + max(insol)^2 *
      t$coefficients[3, ]
    ll <- min(insol) * t$coefficients[2, ] + min(insol)^2 *
      t$coefficients[3, ]
  }


  rmse <- colMeans(t$residuals^2)

  index.neg <- ((ul - ll) <= 0)
  rmse[index.neg] <- 1e+06
  bestfit.index <- which.min(rmse)

  coeff <- t$coefficients[, bestfit.index]


  t.lin <- lm(climate.matrix[, bestfit.index] ~ ins.1)
  coeff.lin <- t.lin$coeff


  rsq1 <- cor(climate.matrix[, bestfit.index], ins.1)^2
  rsq2 <- cor(climate.matrix[, bestfit.index],
              t$fitted.values[, bestfit.index])^2
  coeff.lin <- lm(climate.matrix[, bestfit.index] ~ ins.1)$coeff


  ############# Plotting routine -------
  if (bPlot) {
    at <- basis
    at <- stats::filter(at, rep(1/20, 20), circular = T)
    at.scale <- scale(at)
    insol.scale <- scale(rep(insol, 2))



    at.sc <- attr(at.scale, "scaled:scale")
    at.offset <- attr(at.scale, "scaled:center")

    labels.at <- pretty(c(-3 * at.sc + at.offset, 3 * at.sc + at.offset), 5)
    at.at <- (labels.at - at.offset)/at.sc

    insol.sc <- attr(insol.scale, "scaled:scale")
    insol.offset <- attr(insol.scale, "scaled:center")

    labels.insol <- pretty(c(-3 * insol.sc + insol.offset,
                             3 * insol.sc + insol.offset), 5)
    labels.insol <- labels.insol[labels.insol > 0]
    at.insol <- (labels.insol - insol.offset)/insol.sc

    if (b3plot) {

      par(mfrow = c(2, 2))

      plot(at.scale, ylim = c(-2, 2), axes = F, main = main,
        ylab = clab, xlab = "day of year", type = "l", lwd = 2)
      lines(insol.scale, col = "red", lwd = 2, lty = 2)


      axis(1, at = c(0, 120, 240, 365, 365 + 120, 365 + 240, 365 * 2),
           labels = c(0, 120, 240, 0, 120, 240, 365))
      box()
      axis(2, at = at.at, labels = labels.at)
      axis(4, at = at.insol, labels = labels.insol)


    }
    ############## Plot after cutting ------
    ylab = "insolation (W/m²)"

    plot(at, ylim = c(-2 * at.sc + at.offset, 2 * at.sc + at.offset),
         ylab = clab, xlab = "day of year", type = "l", lwd = 2, axes = F)

    axis(1, at = c(0, 120, 240, 365, 365 + 120, 365 + 240, 365 * 2),
         labels = c(0, 120, 240, 0, 120, 240, 365))
    box()
    axis(2)


    insol.scale <- rep(t.lin$fitted.values, 2)
    iopt <- 365 - 1 * bestfit.index
    insol.scale <- rep(insol.scale[iopt:(iopt + 364)], 2)
    lines(insol.scale, col = "red", lwd = 2, lty = 2)


    insol.scale <- rep(t$fitted.values[, bestfit.index], 2)
    iopt <- 365 - 1 * bestfit.index
    insol.scale <- rep(insol.scale[iopt:(iopt + 364)], 2)
    lines(insol.scale, col = "blue", lwd = 2)


    add <- diff(range(climate.matrix))/10
    main <- ""

    if (b3plot) {
      main <- "response function"
    }

    plot(ins.1, climate.matrix[, bestfit.index], xlab = "insolation (W/m²)",
      ylab = clab, main = main, ylim = range(climate.matrix[, bestfit.index]) +
        c(-add, add))
    lines(ins.1, t$fitted.values[, bestfit.index], col = "blue", lwd = 2)
    lines(ins.1, t.lin$fitted.values, col = "red", lwd = 2)

  }

  return(list(coeff = coeff, coeff.lin = coeff.lin, rsq = c(rsq1,
    rsq2), lag = bestfit.index))

}




# sat.ncep.clim<-read_ncep.clim.day('/Users/tlaepple/data/methaneScaling/data/sat.clim.day.nc',varname='air')
# use_data(sat.ncep.clim)
