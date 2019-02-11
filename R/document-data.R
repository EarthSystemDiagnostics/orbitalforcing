# Document data objects

#' @title Global grid of seasonal cycle amplitudes.
#' @description Global grid of seasonal cycle amplitudes for the past 25 ka calculated using
#' GetSeasCyclePars and transfer functions from EstimateTransferfunctionInsolation usign
#' sat.ncep.clim
#' @format A data frame with 10512 rows and 10 variables:
#' \describe{
#'   \item{\code{lat}}{double COLUMN_DESCRIPTION}
#'   \item{\code{lon}}{double COLUMN_DESCRIPTION}
#'   \item{\code{max.kyear}}{double COLUMN_DESCRIPTION}
#'   \item{\code{max.amp}}{double COLUMN_DESCRIPTION}
#'   \item{\code{min.kyear}}{double COLUMN_DESCRIPTION}
#'   \item{\code{min.amp}}{double COLUMN_DESCRIPTION}
#'   \item{\code{mean.amp}}{double COLUMN_DESCRIPTION}
#'   \item{\code{amp.diff}}{double COLUMN_DESCRIPTION}
#'   \item{\code{sig.sq_c}}{double COLUMN_DESCRIPTION}
#'   \item{\code{sig.sq_a}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
"glob.grid"


#' @title Global grid of insolation to surface temperature transfer functions.
#' @description Global grid of insolation to surface temperature transfer functions.
#' Calculated using EstimateTransferfunctionInsolation and sat.ncep.clim
#' @format A list
#' @details DETAILS
"glob.tf.sat.ncep.clim"

