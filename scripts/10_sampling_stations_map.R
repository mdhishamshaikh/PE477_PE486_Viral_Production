# AIM: To create a map of sampling stationns ####

# 0.0 Setting up #####
source("scripts/0_source.R")

# 1.0 Loading and cropping bathymetry data ####
# Load a 2D NetCDF or GeoTIFF file 
bathymetry <- rast("./Metadata/gebco_2024_n60.0_s50.0_w-2.0_e7.0.nc")  

# Cropping to the North Sea
north_sea_bbox <- ext(-1, 5.5, 52.5, 56)  
bathymetry_north_sea <- crop(bathymetry, north_sea_bbox)

# Converting to data frame for ggplot2
bathymetry_df <- as.data.frame(bathymetry_north_sea, xy = TRUE)
colnames(bathymetry_df) <- c("Longitude", "Latitude", "Depth")

# Updating the bathymetry data
bathymetry_df$Depth <- ifelse(bathymetry_df$Depth > 0, NA, bathymetry_df$Depth)  # Setting land (above 0) to NA
bathymetry_df$Depth <- -bathymetry_df$Depth  # Reversing depth (more negative = more depth)


# 2.0 Importing data for PE crusies ####

pe_df <- read.csv( "./results/PE477_PE486_3depths_combined.csv") %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_")) # No VP assay was performed

# Extracting stations and their coordinates

stations <- pe_df %>%
  dplyr::select(c(Location, Station_Number, Latitude, Longitude, Season)) %>%
  drop_na() %>%
  dplyr::mutate(Station_Number = round(Station_Number)) %>%
  distinct() 

# This gives 12 stations 
# Remember we have 12.1 an 12.2 that got rounded off to 12. These samples were taken a day apart.

# Assigning colors and shapes ti seasons 

custom_season_colors <- c(`Autumn (Sept 2020)` = "#E69F00",
                          `Spring (Apr 2021)` = "#009E73")
custom_season_fills <- c(`Autumn (Sept 2020)` = "#003049",
                          `Spring (Apr 2021)` = "#DC2828")
custom_season_shapes <- c(`Autumn (Sept 2020)` = 19, 
                          `Spring (Apr 2021)` = 17)

# 3.0 Plotting stations map####
stations_plot<- ggplot() +
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
    na.value = "black",  # Land areas in black
    guide = guide_colorbar(reverse = TRUE) 
  ) +
  geom_point(data = stations, aes(x = Longitude, y = Latitude, shape = Season#, color = Season
                                  ),
             size = 5, alpha = 0.8,  
             ) +  # Plot VBR points
  #geom_point(data = stations, aes(x = Longitude, y = Latitude))+
  geom_text(data = stations, aes(x = Longitude, y = Latitude, label = Station_Number), 
            color = "black", size = 5, nudge_y = -0.15, nudge_x = -0.1, check_overlap = F) + 
  # geom_text_repel(data = stations, aes(label = Station_Number), size = 3, box.padding = 0.3, max.overlaps = 100) +
  #scale_colour_gradientn(colours = (ODV_colours)) +
  # coord_cartesian(expand = F)+
 # facet_grid(~ Location) +
  scale_x_continuous(
    #name = "Longitude",
    labels = function(x) {
      ifelse(x < 0, paste0(abs(x), " °W"),
             ifelse(x > 0, paste0(x, " °E"), "0 °"))
    },
    position = "top") +
  scale_y_continuous(
    #name = "Latitude (°N)",
    labels = function(y) paste0(abs(y), " °N"),
    position = "right"
  ) +
  #scale_color_manual(values = custom_season_colors) +
  scale_shape_manual(values = custom_season_shapes) +
  labs(
   # title = "Bathymetry of the North Sea",
    x = NULL,
    y = NULL,
    shape = "Seasons"
  ) +
  # annotation_scale(
  #   location = "bottomleft", 
  #   width_hint = 0.2, 
  #   height = unit(0.3, "cm"),
  #   text_col = "black"
  # ) + 
  # annotation_north_arrow(
  #   location = "topright",
  #   which_north = "true",
  #   pad_x = unit(0.5, "cm"),
  #   pad_y = unit(0.5, "cm"),
  #   style = north_arrow_fancy_orienteering()
  # ) +
  coord_fixed(ratio = 1.3,
              expand = F
  ) +  # Preserve aspect ratio
  theme_bw(base_size = 15) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.border = element_rect(linewidth = 2, color = "black"),
    plot.margin = margin(10, 10, 10, 20)
  )
stations_plot
ggsave(stations_plot, filename = "./figures/stations_plot_with_bathymetry.svg", width = 10, height =8, dpi = 800)




# Define the bounding box for the North Sea bathymetry map
bbox_europe <- data.frame(
  xmin = -1, xmax = 5.5, ymin = 52.5, ymax = 56
)

# Get a simple world map
map_data <- map_data("world")

# Plot the map with bounding box
map_plot <- ggplot() +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), fill = "grey50", color = "black") +
  geom_rect(data = bbox_europe, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", size = 1.2) +  # Bounding box
  coord_fixed(xlim = c(-10, 10), ylim = c(35, 65), expand = FALSE, ratio = 1.3) +
  scale_x_continuous(
    #name = "Longitude",
    labels = function(x) {
      ifelse(x < 0, paste0(abs(x), " °W"),
             ifelse(x > 0, paste0(x, " °E"), "0 °"))
    }) +
  scale_y_continuous(
    #name = "Latitude (°N)",
    labels = function(y) paste0(abs(y), " °N")
  ) +
  theme_minimal() +
  labs(#title = "Location of North Sea Bathymetry Study Area",
       #subtitle = "Bounding box shows limits of detailed bathymetry map",
       x = NULL,
       y = NULL)+
 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        panel.border = element_rect( color = "black"))
map_plot
ggsave(map_plot, filename = "./figures/europe_map.svg", width = 10, height =8, dpi = 800)

