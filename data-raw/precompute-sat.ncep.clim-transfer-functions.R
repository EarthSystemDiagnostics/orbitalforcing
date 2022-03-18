## Get transfer functions from sat.ncep.clim
library(orbitalforcing)
library(tidyverse)

clim <- sat.ncep.clim

GetTF <- function(df){
  climate.modern <- pfields::SelSpace3D(clim,
                               lat1 = df$lat, lon1 = df$lon) - 273.15

tf <- orbitalforcing::EstimateTransferfunctionInsolation(
  climate = climate.modern, latitude = df$lat)

return(list(tf))

}


glob.tf.sat.ncep.clim <-
  crossing(lat = attr(clim, "lat"),
            lon = attr(clim, "lon")) %>%
  group_by(lat, lon) %>%
  plyr::dlply(., c("lat", "lon"), GetTF, .progress = "text")

lon1 <- attr(glob.tf.sat.ncep.clim, "split_labels")$lon
attr(glob.tf.sat.ncep.clim, "split_labels")$lon <- ifelse(lon1 >=180, lon1 - 360, lon1)




## precompute global grid of amps

clim <- sat.ncep.clim

glob.grid <-
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



#Load the matrix contains data from Berger and Loutre (1991),
##'   downloaded as ORBIT91 from ncdc.noaa.gov
#orbit91<-read.table(paste(path,"ins_data.txt",sep=""))
#orbital_global<-orbital_parameters((0:50000)/10)


usethis::use_data(orbital_global, overwrite = TRUE)
usethis::use_data(orbit91, overwrite = T)
usethis::use_data(glob.tf.sat.ncep.clim, overwrite = T)
usethis::use_data(glob.grid, overwrite = T)
usethis::use_data(sat.ncep.clim, overwrite = T)

usethis::use_data(#glob.grid, orbit91,
  orbital_global,
  internal = TRUE, overwrite = TRUE)


# load("data/orbital_global.rda")
# orb_g_data <- orbital_global
#
# rm(orbital_global)
#
# load("R/sysdata.rda")
# orb_g_int <- orbital_global
#
# lapply(1:3, function(x) {
#   table(orb_g_data[[x]] == orb_g_int[[x]])
#})
