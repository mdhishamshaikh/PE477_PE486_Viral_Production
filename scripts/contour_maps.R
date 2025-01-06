library(tidyverse)

ODV_colours <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")


url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
str(abundance)

#Extracting first entries per Location/Station_Number combinations, as some stations had multiple depths and this information is missing from the fiel on GitHub.
#The first one is the depth for VP assays.
abundance<- abundance %>%
  dplyr::filter(Location %in% c("PE477", "PE486")) %>%
  dplyr::select(-c(Temperature, Salinity, ends_with("Sample_Name"), Expt_Date, TON))


# Example: Switch to Temperature if Total_Bacteria has no variation
ggplot(data_7m, aes(x = Longitude, y = Latitude, z = VBR)) +
  #geom_contour_filled(aes(fill = after_stat(level)), bins = 10) +
  labs(
    title = "Contour Map of Temperature at 7 m Depth",
    x = "Longitude",
    y = "Latitude",
    fill = "Temperature"
  ) +
  theme_minimal(base_size = 14)
# install.packages("ggplot2")
library(ggplot2)

ggplot(abundance, aes(x = Longitude, y = Latitude, fill = Total_Bacteria)) +
  geom_tile() +  # Create the grid
  scale_y_reverse() +  # Reverse depth axis to mimic ODV-style
  scale_fill_viridis_c() +  # Add a perceptually uniform color scale
  labs(
    title = "Bacterial Abundance Across Depth and Station",
    x = "Station Number",
    y = "Depth (m)",
    fill = "Bacterial Abundance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

##################

#https://theoceancode.netlify.app/post/odv_bathy/
  library(tidyverse)
library(lubridate)
library(reshape2)
library(MBA)
library(mgcv)
library(marmap)
library(FNN)

data_7m <- abundance %>%
  dplyr::filter(Depth == 7)

ggplot(data = data_7m, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(colour = VBR), size = 20, alpha = 0.75) +
  geom_point()+
  facet_grid(. ~ Location) +
  scale_colour_gradientn(colours = rev(ODV_colours))


library(dplyr)
library(MBA)
library(reshape2)
library(ggplot2)

# Define the global extent of coordinates
global_extent <- list(
  Longitude = seq(min(data_7m$Longitude), max(data_7m$Longitude), length.out = 300),
  Latitude = seq(min(data_7m$Latitude), max(data_7m$Latitude), length.out = 300)
)

# Function to interpolate for a single location with the same extent
interpolate_with_common_extent <- function(location_data, extent) {
  # Perform MBA interpolation
  mba_result <- mba.surf(
    location_data[c("Longitude", "Latitude", "VBR")],
    no.X = length(extent$Longitude),
    no.Y = length(extent$Latitude),
    extend = TRUE
  )
  
  # Set the common extent to the interpolation result
  dimnames(mba_result$xyz.est$z) <- list(extent$Longitude, extent$Latitude)
  melted_data <- melt(mba_result$xyz.est$z, varnames = c('Longitude', 'Latitude'), value.name = 'VBR') %>%
    mutate(
      VBR = round(VBR, 5),
      Longitude = as.numeric(as.character(Longitude)),
      Latitude = as.numeric(as.character(Latitude))
    )
  
  # Add the Location back
  melted_data$Location <- unique(location_data$Location)
  return(melted_data)
}

# Apply the function for each location and combine results
ctd_mba <- data_7m %>%
  group_split(Location) %>%
  lapply(interpolate_with_common_extent, extent = global_extent) %>%
  bind_rows()

library(ggplot2)
library(maps)

# Define limits from your data
lon_range <- range(data_7m$Longitude, na.rm = TRUE)
lat_range <- range(data_7m$Latitude, na.rm = TRUE)

# Plot with country borders and explicit limits
ggplot(data = ctd_mba, aes(x = Longitude, y = Latitude)) +
  geom_raster(aes(fill = VBR)) +
  geom_point(data = data_7m, aes(x = Longitude, y = Latitude), colour = "black", size = 1) +
  geom_contour(aes(z = VBR), binwidth = 2, colour = "black", alpha = 0.2) +
  borders("world", colour = "black", size = 0.5) +  # Add country borders
  facet_wrap(~ Location) +  # Facet by Location
  scale_fill_gradientn(colours = rev(ODV_colours)) +
  labs(
    y = "Latitude",
    x = "Longitude",
    fill = "VBR"
  ) +
  coord_cartesian(
    xlim = lon_range,  # Set Longitude limits
    ylim = lat_range   # Set Latitude limits
  ) +
  theme_minimal()

##########


library(dplyr)
library(MBA)
library(reshape2)
library(ggplot2)

# Define the global extent of coordinates with 5-degree extensions
global_extent <- list(
  Longitude = seq(min(data_7m$Longitude, na.rm = TRUE) - 1, max(data_7m$Longitude, na.rm = TRUE) + 1, length.out = 300),
  Latitude = seq(min(data_7m$Latitude, na.rm = TRUE) - 1, max(data_7m$Latitude, na.rm = TRUE) + 1, length.out = 300)
)

# Function to interpolate for a single location with the extended extent
interpolate_with_common_extent <- function(location_data, extent) {
  # Perform MBA interpolation with extended limits
  mba_result <- mba.surf(
    location_data[c("Longitude", "Latitude", "VBR")],
    no.X = length(extent$Longitude),
    no.Y = length(extent$Latitude),
    extend = TRUE
  )
  
  # Set the extended extent to the interpolation result
  dimnames(mba_result$xyz.est$z) <- list(extent$Longitude, extent$Latitude)
  melted_data <- melt(mba_result$xyz.est$z, varnames = c('Longitude', 'Latitude'), value.name = 'VBR') %>%
    mutate(
      VBR = round(VBR, 5),
      Longitude = as.numeric(as.character(Longitude)),
      Latitude = as.numeric(as.character(Latitude))
    )
  
  # Add the Location back
  melted_data$Location <- unique(location_data$Location)
  return(melted_data)
}

# Apply the function for each location and combine results
ctd_mba <- data_7m %>%
  group_split(Location) %>%
  lapply(interpolate_with_common_extent, extent = global_extent) %>%
  bind_rows()

# Plot with the extended interpolation results
ggplot(data = ctd_mba, aes(x = Longitude, y = Latitude)) +
  geom_raster(aes(fill = VBR)) +
  geom_point(data = data_7m, aes(x = Longitude, y = Latitude), colour = "black", size = 1) +
  geom_contour(aes(z = VBR), binwidth = 2, colour = "black", alpha = 0.2) +
  borders("world", colour = "black", fill = "black", size = 0.5) +  # Add country borders
  facet_wrap(~ Location) +  # Facet by Location
  scale_fill_gradientn(colours = rev(ODV_colours)) +
  labs(
    y = "Latitude",
    x = "Longitude",
    fill = "VBR"
  ) +
  coord_cartesian(
    xlim = range(global_extent$Longitude),  # Use extended Longitude limits
    ylim = range(global_extent$Latitude)   # Use extended Latitude limits
  ) +
  theme_test()


# With OCE #####

# Load necessary packages
library(oce)
library(akima)


# Interpolate to create a grid
interp_result <- akima::interp(
  x = data_7m$Longitude,
  y = data_7m$Latitude,
  z = data_7m$Nitrate,
  xo = seq(min(data_7m$Longitude), max(data_7m$Longitude), length = 100),
  yo = seq(min(data_7m$Latitude), max(data_7m$Latitude), length = 100)
)

# Plotting in ODV style
oce::imagep(data_7m$Longitude, data_7m$Latitude, data_7m$Nitrate,
            col = oceColorsJet, 
            zlab = "Nitrate (µM)",
            xlab = "Longitude", 
            ylab = "Latitude",
            main = "Nitrate Distribution")

# Add station points
points(data$longitude, data$latitude, pch = 21, bg = "white", cex = 1.5)

global_extent <- list(
  Longitude = seq(min(data_7m$Longitude), max(data_7m$Longitude), length.out = 300),
  Latitude = seq(min(data_7m$Latitude), max(data_7m$Latitude), length.out = 300)
)

# Geom point plots on map
ggplot(data = data_7m, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(colour = VBR), size = 20, alpha = 0.75) +
  geom_point()+
  facet_grid(. ~ Location) +
  scale_colour_gradientn(colours = rev(ODV_colours)) +
  borders("world", colour = "black", size = 0.5) +  # Add country borders
  labs(
    y = "Latitude",
    x = "Longitude",
    fill = "VBR"
  ) +
  coord_cartesian(
    xlim = range(global_extent$Longitude),  # Use extended Longitude limits
    ylim = range(global_extent$Latitude)   # Use extended Latitude limits
  ) +
  theme_test()




# Define global_extent for limits
global_extent <- data.frame(
  Longitude = c(-10, 10),  # Adjust these limits as needed
  Latitude = c(50, 60)     # Adjust these limits as needed
)

# Custom ODV-like color palette (replace with your actual colors)
ODV_colours <- oce::oceColorsJet(100)

# Plot with black land fill
ggplot(data = data_7m, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(colour = VBR), size = 20, alpha = 0.75) +
  geom_point()+scale_colour_gradientn(colours = rev(ODV_colours)) +     # Reverse color palette
  borders("world", fill = "black", colour = "black") +     # Add black land
  labs(
    y = "Latitude",
    x = "Longitude",
    colour = "VBR"
  ) +
  coord_cartesian(
    xlim = range(global_extent$Longitude),  # Set Longitude limits
    ylim = range(global_extent$Latitude)   # Set Latitude limits
  ) +
  facet_grid(. ~ Location) +  # Facet by Location
  theme_test() +
  theme(
    panel.background = element_rect(fill = "white"),  # Ocean background color
    panel.grid = element_line(colour = "white")           # Optional grid lines
  )


# Plot with Earth's curvature
ggplot() +
  geom_sf(data = world, fill = "black", color = "black") +  # Land in black
  geom_point(data = data_7m, aes(x = Longitude, y = Latitude, colour = VBR),
             size = 5, alpha = 0.75) +  # Plot VBR points
  scale_colour_gradientn(colours = rev(ODV_colours)) +  # Reverse color palette
  labs(
    y = "Latitude",
    x = "Longitude",
    colour = "VBR"
  ) +
  coord_sf(crs = "+proj=moll",  # Mollweide projection
           xlim = range(data_7m$Longitude) + c(-2, 2),  # Extend Longitude limits
           ylim = range(data_7m$Latitude) + c(-5, 5)) + # Extend Latitude limits
  facet_grid(. ~ Location) +  # Facet by Location
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),  # Ocean background color
    panel.grid = element_line(colour = "gray")        # Optional grid lines
  )



# High resolution North Sea bathymetry plot


library(terra)
library(ggplot2)

# Load a 2D NetCDF or GeoTIFF file (replace with your file path)
bathymetry <- rast("./Metadata/gebco_2024_n60.0_s50.0_w-2.0_e7.0.nc")  # or .tif

# Crop to your region of interest (e.g., North Sea)
north_sea_bbox <- ext(-1, 5.5, 52.5, 56)  # Adjust limits as needed
bathymetry_north_sea <- crop(bathymetry, north_sea_bbox)

# Convert to data frame for ggplot2
bathymetry_df <- as.data.frame(bathymetry_north_sea, xy = TRUE)
colnames(bathymetry_df) <- c("Longitude", "Latitude", "Depth")

# Update the bathymetry data
bathymetry_df$Depth <- ifelse(bathymetry_df$Depth > 0, NA, bathymetry_df$Depth)  # Set land (above 0) to NA
bathymetry_df$Depth <- -bathymetry_df$Depth  # Reverse depth (more negative = more depth)

# Plot bathymetry with land in black
ggplot() +
  geom_raster(data = bathymetry_df, aes(x = Longitude, y = Latitude, fill = Depth), na.rm = TRUE) +
  scale_fill_gradientn(
    colours = rev(c(
      "white",       # Land (0 m)
      "#E3F2FD",     # Lightest blue (-5 m)
      #"#BBDEFB",     # Very light blue (-10 m)
      "#90CAF9",     # Light blue (-15 m)
      #"#64B5F6",     # Light sky blue (-20 m)
      "#42A5F5",     # Light dodger blue (-25 m)
      #"#2196F3",     # Dodger blue (-30 m)
      # "#1E88E5",     # Medium blue (-35 m)
      # "#1976D2",     # Steel blue (-40 m)
      # "#1565C0",     # Royal blue (-50 m)
      "#0D47A1",     # Deeper royal blue (-60 m)
      "#0B408A",     # Darker blue (-70 m)
      # "#093974",     # Navy blue (-80 m)
      "#082E63",     # Deep navy blue (-90 m)
      # "#061D4F",     # Indigo blue (-100 m)
      "#04143D",     # Very dark indigo (-125 m)
      # "#030C2C",     # Almost black blue (-150 m)
      # "#020822",     # Very deep blue (-175 m)
      "#010417",     # Deepest blue (-200 m)
      "black"        # Very dark blue/black (-250 m)
    )),
    values = scales::rescale(c(0, -5, -10, -15, -20, -25, -30, 
                               -35, -40, -45, -50
                               # , -55, -60,
                               # -70, -75, -80, -100, -130, -170, -250
                               )),  # Fine-grained depth breaks
    name = "Depth (m)", 
    na.value = "black"  # Land areas in black
  ) +
  geom_point(data = data_7m, aes(x = Longitude, y = Latitude, colour = VBR),
             size = 30, alpha = 0.8, shape = 19) +  # Plot VBR points
  geom_point(data = data_7m, aes(x = Longitude, y = Latitude))+
  geom_text(data = data_7m, aes(x = Longitude, y = Latitude, label = Station_Number), 
            color = "black", size = 5, nudge_y = -0.1, nudge_x = 0.1, check_overlap = F) + 
  scale_colour_gradientn(colours = (ODV_colours)) +
# coord_cartesian(expand = F)+
  facet_grid(~ Location) +
  labs(
    title = "Bathymetry of the North Sea",
    x = "Longitude",
    y = "Latitude"
  ) +
  coord_fixed(ratio = 1.3,
              expand = F
              ) +  # Preserve aspect ratio
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


#Lambert projection



# Saving 


map_data <- read.csv("./results/ctd_profiles_abundance_nutrients_3_depths.csv")
