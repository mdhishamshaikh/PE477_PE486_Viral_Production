# Lets plot Temeprature on a bathymetry m,map

library(marmap)
library(ggplot2)


t_df <- vp_pe_df %>%
  select(Latitude, Longitude, Temperature)


# Defining boundaries
lon_min <- min(t_df$Longitude) - 10
lon_max <- max(t_df$Longitude) + 10
lat_min <- min(t_df$Latitude) - 10
lat_max <- max(t_df$Latitude) + 10

#Fetchin g bathymetry data from NOAA
bathy_data <- getNOAA.bathy(lon1 = lon_min, lon2 = lon_max, lat1 = lat_min, lat2 = lat_max, resolution = 1)

plot(bathy_data, image = TRUE, deep = -100, shallow = 0, step = 10, land = TRUE, n = 100)

with(t_df, {
  points(Longitude, Latitude, pch = 20, 
         col = ifelse(is.na(Temperature), "grey", 
                      colorRampPalette(c("blue", "red"))((Temperature - min(Temperature, na.rm = TRUE)) / 
                                                           (max(Temperature, na.rm = TRUE) - min(Temperature, na.rm = TRUE)) * length(Temperature))), 
         cex = 1.5)
  text(Longitude, Latitude, labels = ifelse(is.na(Temperature), "", round(Temperature, 1)), cex = 0.8, pos = 3)
})
