## Get transfer functions from sat.ncep.clim
library(orbitalforcing)
library(tidyverse)

clim <- sat.ncep.clim

GetTF <- function(df){
  climate.modern <- ecustools::SelSpace3D(clim,
                                          lat1 = df$lat, lon1 = df$lon) - 273.15

tf <- orbitalforcing::EstimateTransferfunctionInsolation(
  climate = climate.modern, latitude = df$lat)

return(list(tf))

}


glob.tf.sat.ncep.clim <-
  crossing(lat = attr(clim, "lat"),
            lon = attr(clim, "lon")) %>%
  #crossing(lat = c(23, 56),
   #        lon = c(78, 89)) %>%
  group_by(lat, lon) %>%
  plyr::dlply(., c("lat", "lon"), GetTF, .progress = "text")

lon1 <- attr(glob.tf.sat.ncep.clim, "split_labels")$lon
attr(glob.tf.sat.ncep.clim, "split_labels")$lon <- ifelse(lon1 >=180, lon1 - 360, lon1)

devtools::use_data(glob.tf.sat.ncep.clim, overwrite = T)



## precompute global grid of amps

clim <- sat.ncep.clim

glob.grid <-
  # crossing(lat = attr(clim, "lat"),
  #          lon = ifelse(attr(clim, "lon") >= 180,
  #                       attr(clim, "lon") -360, attr(clim, "lon"))) %>%
   crossing(lat = attr(clim, "lat"),
            lon = attr(clim, "lon")) %>%
  group_by(lat, lon) %>%
  do({
    pars <- GetSeasCyclePars(lon = .$lon, lat=.$lat)
    as_tibble(pars)
  })


glob.grid <- glob.grid %>%
  ungroup() %>%
  mutate(lon = ifelse(lon >= 180,
                lon -360, lon))

devtools::use_data(glob.grid, overwrite = T)

