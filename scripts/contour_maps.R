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
# 0.0 Set up ####

# Loading necessary libraries 
library(tidyverse)
library(terra)
library(cowplot)
library(colorspace)




# Loading the csv file with CTD, nutrients, and abundance data from 7 m, 15 m, and 30 m.
map_data <- read.csv("./results/ctd_profiles_abundance_nutrients_3_depths.csv")


# Defining custom color palette
custom_color_palette <- c(
  "PE477_1" = "#1f77b4", "PE477_2" = "#ff7f0e", "PE477_3" = "#2ca02c",
  "PE477_4" = "#d62728", "PE477_5" = "#9467bd", "PE477_6" = "#8c564b",
  "PE486_1" = "#e377c2", "PE486_2" = "#7f7f7f", "PE486_3" = "#bcbd22",
  "PE486_4" = "#17becf", "PE486_5" = "#aec7e8", "PE486_6" = "#ffbb78",
  "PE477_7" = "#98df8a", "PE486_7" = "#ff9896",
  "PE477" = colorspace::lighten("#0c1844", 0.0), "PE486" = colorspace::lighten("#850000", 0.2)
)
ODV_colours <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")



# 1.0 Bathymetry map of the North Sea with sampling stations ####

# Loading a 2D NetCDF downloaded from GEBCO
bathymetry <- rast("./Metadata/gebco_2024_n60.0_s50.0_w-2.0_e7.0.nc") 

# Cropping to the North Sea
north_sea_bbox <- ext(-1.5, 6, 52, 56.5)  
bathymetry_north_sea <- crop(bathymetry, north_sea_bbox)

# Converting to data frame for ggplot2
bathymetry_df <- as.data.frame(bathymetry_north_sea, xy = TRUE)
colnames(bathymetry_df) <- c("Longitude", "Latitude", "Depth")

# Setting all above sea surface (> 0 m) values to NA (land)
bathymetry_df$Depth <- ifelse(bathymetry_df$Depth > 0, NA, bathymetry_df$Depth)  

# Reversing depth
bathymetry_df$Depth <- -bathymetry_df$Depth 


# Subsettig 7 m data for ease

data_7m <- map_data %>%
  dplyr::filter(Depth == 7)

# Plotting bathymetry
North_Sea_bathymetry_sampling_stations_plot<-  ggplot() +
    geom_raster(data = bathymetry_df, aes(x = Longitude, y = Latitude, fill = Depth), na.rm = TRUE) +
    scale_fill_gradientn(
      colours = rev(c(
        "white",       
        "#E3F2FD",    
        "#90CAF9",
        "#42A5F5",    
        "#0D47A1",   
        "#0B408A",   
        "#082E63",    
        "#04143D",     
        "#020822",     
        "#010417",    
        "black"       
      )),
      values = scales::rescale(c(0, -5, -10, -15, -20, -25, 
                                 -30, -35, -40, -45, -50)), 
      name = "Depth (m)", 
      na.value = "black",
      guide = guide_colorbar(reverse = TRUE)
    ) +
   geom_point(data = data_7m, aes(x = Longitude, y = Latitude, color = Location, shape = Location),
              size = 3)+
    geom_text(data = data_7m, aes(x = Longitude, y = Latitude, label = Station_Number), 
              color = "black", size = 5, nudge_y = 0.15, nudge_x = 0.1, check_overlap = F) + 
    +scale_colour_gradientn(colours = rev(ODV_colours))  +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    coord_fixed(ratio = 1,
                expand = F
    ) + 
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

North_Sea_bathymetry_sampling_stations_plot

ggsave("./figures/North_Sea_bathymetry_sampling_stations.svg", plot = North_Sea_bathymetry_sampling_stations_plot,
       dpi = 800, width = 8.27,  units = "in")


# 2.0 Plotting biotic and abiotic variables on North sea bathymetry ####


plot_bathymetry_cowplot <- function(bathymetry_df, data, depths, custom_palette, variable_groups, output_prefix = "bathymetry_plots") {
  
  # Looping through each depth
  for (depth in depths) {
   
    data_depth <- data %>%
      filter(Depth == depth)
    
    # Looping through variable groups 
    for (i in seq_along(variable_groups)) {
      variables_to_plot <- variable_groups[[i]]
      
      # Generating individual plots for each variable
      plot_list <- lapply(variables_to_plot, function(var) {
        # Pivot longer for the current variable
        # data_long <- data_depth %>%
        #   select(Location_Station_Number, Location, Station_Number, Latitude, Longitude, Depth, all_of(var)) %>%
        #   pivot_longer(cols = all_of(var), names_to = "variable", values_to = "value")
        # 
        # Generate the plot for the current variable
        ggplot() +
          # Raster for bathymetry
          geom_raster(data = bathymetry_df, aes(x = Longitude, y = Latitude, fill = Depth), na.rm = TRUE) +
          # Gradient fill for Depth
          scale_fill_gradientn(
            colours = rev(c(
              "white", "#E3F2FD", "#90CAF9", "#42A5F5", "#0D47A1", 
              "#0B408A", "#082E63", "#04143D", "#020822", "#010417", "black"
            )),
            values = scales::rescale(c(0, -5, -10, -15, -20, -25, 
                                       -30, -35, -40, -45, -50)),
            name = "Depth (m)",
            na.value = "black",
            guide = guide_colorbar(reverse = TRUE)
          ) +
          guides(fill = "none") +
          # Points for data
          geom_point(data = data_depth, aes(x = Longitude, y = Latitude, colour = !!sym(var)),
                     size = 10, alpha = 0.8, shape = 19) +
          # Station numbers as text labels
          geom_text(data = data_depth, aes(x = Longitude, y = Latitude, label = Station_Number), 
                    color = "black", size = 4, nudge_y = 0, nudge_x = 0, check_overlap = FALSE) + 
          # Color scale for variable value
          scale_colour_gradientn(colours = rev(custom_palette)) +
          # Facet by Location
          facet_wrap(. ~ Location, ncol = 2, nrow = 2) +
          labs(
            x = "Longitude",
            y = "Latitude",
            title = paste("Depth:", depth, "m |", var)
          ) +
          coord_fixed(ratio = 1, expand = FALSE) + 
          theme_bw(base_size = 14) +
          theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white")
          )
      })
      
      # Combine the plots for the current variable group
      combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)

      # Save the combined plot as an SVG with dynamic width and height
      n_plots <- length(plot_list)  # Number of plots in the current group
      n_cols <- 2  # Number of columns in the grid
      n_rows <- ceiling(n_plots / n_cols)  # Calculate the number of rows based on plots and columns

      ggsave(
        filename = paste0("./figures/", output_prefix, "_depth_", depth, "_group_", i, ".svg"),
        plot = combined_plot,
        width = 6 * n_cols,  # Adjust width (e.g., 6 inches per column)
        height = 2.5 * n_rows,  # Adjust height (e.g., 5 inches per row)
        dpi = 300
      )
      # # Save the combined plot as an SVG
      # ggsave(
      #   filename = paste0(output_prefix, "_depth_", depth, "_group_", i, ".svg"),
      #   plot = combined_plot,
      #   width = 12, height = 10, dpi = 300
      # )
    }
  }
}

# To create split plots, adjust this list 
variable_groups <- list(
  c( "Temperature", "Salinity", "Density", "Conductivity", "Turbidity", "Oxygen", "Fluorescence",
      "Nitrate", "Nitrite", "Silicate", "Phosphate", "VBR", "Total_Bacteria", "Total_Viruses")
  #c("Total_Bacteria", "Total_Viruses", "VBR", "Nitrate", "Nitrite", "Silicate", "Phosphate")
)

# Depths to plot
depths <- unique(map_data$Depth)

# Adjusting Total Bacteria and Total Viruses to (in millions)
map_data_adj_counts  <- map_data %>%
  mutate(across(c(Total_Bacteria, HNA, LNA, Total_Viruses, V1, V2, V3), ~ .x / 1e+6))

# Plotting
plot_bathymetry_cowplot(
  bathymetry_df = bathymetry_df,
  data = map_data_adj_counts,
  depths = depths,
  custom_palette = ODV_colours,
  variable_groups = variable_groups,
  output_prefix = "bathymetry_plot"
)



