# AIM: To extract CTD profiles for all PE477 & PE486 stations, and extract values for the sampling depths (7, 15, and 30 m)

# 0.0 Setting up ####
source("./scripts/0_source.R")

# 1.0 Importing CTD profiles and selecting stations ####
file_path <- "./data/ctd_profiles"
file_list <- list.files(path = file_path, pattern = "*.csv", full.names = TRUE)

# Function to add file name as a column
read_and_label <- function(file) {
  
  df <- read.csv(file)
  df$file_name <- basename(file)
  return(df)
}


# Changing file names and keeping only the stations needed.
combined_df <- bind_rows(lapply(file_list, read_and_label))
combined_df$file_name <- sub("\\.csv$", "", combined_df$file_name)
combined_df <- combined_df[combined_df$file_name %in% c("PE486_S02C01", "PE486_S04C02", "PE486_S06C01", "PE486_S08C01", 
                                                        "PE486_S09C01", "PE486_S10C01", "PE486_S11C01", "PE486_S12C01", 
                                                        "64PE477_S00C01", "64PE477_S04C01", "64PE477_S05C01", "64PE477_S07C01", 
                                                        "64PE477_S10C01", "64PE477_S15C01", "64PE477_S17C01"), ]

# Simplifying station names and adding Location, Station_Number columns
name_mapping <- c(
  "PE486_S02C01" = "PE486_1",
  "PE486_S04C02" = "PE486_2",
  "PE486_S06C01" = "PE486_3",
  "PE486_S08C01" = "PE486_4",
  "PE486_S09C01" = "PE486_5",
  "PE486_S10C01" = "PE486_6",
  "PE486_S11C01" = "PE486_7",
  "PE486_S12C01" = "PE486_8",
  "64PE477_S00C01" = "PE477_1",
  "64PE477_S04C01" = "PE477_2",
  "64PE477_S05C01" = "PE477_3",
  "64PE477_S07C01" = "PE477_4",
  "64PE477_S10C01" = "PE477_5",
  "64PE477_S15C01" = "PE477_6",
  "64PE477_S17C01" = "PE477_7"
)

combined_df$Location_Station_Number <- name_mapping[combined_df$file_name]
combined_df <- combined_df %>%
  separate(Location_Station_Number, into = c("Location", "Station_Number"), sep = "_", remove = F)

unique(combined_df$Location_Station_Number)


# Removing station PE477_7 and PE486_8 as they do not have viral production data attached to them. 
# Will need these stations for 'omics work though.
# combined_df <- combined_df %>%
#   dplyr::filter(!Location_Station_Number %in% c("PE477_7", "PE486_8")) 


# 2.0 Selecting variables and simplifying names ####

combined_df <- combined_df %>%
  rename(
    Original_station_cast = file_name,
    Depth = depSM,
    Salinity = sal00,
    Temperature = t090C,
    Pressure = prDM,
    Conductivity = C0,
    Oxygen = sbeox0Mm.L,
    Turbidity = turb,
    Chlorophyll = flC
  )

# 3.0 Calculating density and removing surface values values ####

# Gibbs SeaWater 
combined_df$Density <- gsw::gsw_rho(combined_df$Salinity, combined_df$Temperature, combined_df$Pressure)

combined_df <- combined_df %>%
  dplyr::filter(Depth > 3) # Making sure first 3 meters are not considered as the disturbance from dropping the CTD frame creates large fluctuations in measurements.

ctd_profiles<- combined_df %>%
  dplyr::select(c("Original_station_cast", "Location_Station_Number", "Depth", "Salinity", "Temperature", "Pressure", "Density", "Conductivity", "Oxygen", "Turbidity", "Chlorophyll")) %>%
  separate(Location_Station_Number, c("Location", "Station_Number"), "_", remove = F)


# Max depth
ctd_profiles <- ctd_profiles %>%
  group_by(Location_Station_Number) %>%
  mutate(Max_Depth = max(Depth)) %>%
  ungroup()

# Writing output as a csv  
write.csv(ctd_profiles, "./results/ctd_profiles/ctd_profiles.csv", row.names = F)

# 4.0 Reshaping and plotting #####
# Reading CTD profiles csv
ctd_profiles <- read.csv("./results/ctd_profiles/ctd_profiles.csv")

ctd_long_df <- ctd_profiles %>%
  dplyr::select(-Location, -Station_Number) %>%  
  reshape2::melt(id.vars = c("Original_station_cast", "Location_Station_Number", "Depth")) %>%  
  mutate(variable = factor(variable, levels = c(
    "Salinity", "Temperature", "Pressure", "Density", 
    "Conductivity", "Oxygen", "Turbidity", "Chlorophyll"
  )))


unique(ctd_long_df$variable)

# Plotting
pdf("./results/ctd_profiles/depth_profiles_by_station.pdf", width = 10, height = 7)


for (station in unique(ctd_long_df$Location_Station_Number)) {
  
  
  station_df <- subset(ctd_long_df, Location_Station_Number == station)
  
  p <- ggplot(station_df, aes(x = value, y = Depth)) +
    geom_point(size = 0.1) +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) + 
    labs(title = paste("Depth Profiles for", station),
         x = "Parameter Value",
         y = "Depth (m)") +
    scale_y_reverse() +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          strip.text = element_text(size = 8))              
  
  print(p)
}

dev.off()

# As images

# Creating a folder for the images
dir.create("./results/ctd_profiles/per_station", showWarnings = FALSE)

for (station in unique(ctd_long_df$Location_Station_Number)) {
  
  
  station_df <- subset(ctd_long_df, Location_Station_Number == station)
  
  p <- ggplot(station_df, aes(x = value, y = Depth)) +
    geom_point(size = 0.1) +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +  
    labs(title = paste("Depth Profiles for", station),
         x = "Parameter Value",
         y = "Depth (m)") +
    scale_y_reverse() +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          strip.text = element_text(size = 8))             
  
 
  ggsave(
    filename = paste0("./results/ctd_profiles/per_station/depth_profile_", station, ".png"),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300  
  )
}


# 5.0 Extracing CTD variables for 7, 15 and 30 m ####

ctd_profiles <- read.csv("./results/ctd_profiles/ctd_profiles.csv")

# Target depths
target_depths <- c(7, 15, 30) # Some station might not have 30 m. Use the distance variable to figure out. Maybe worth adding max depth

# Function to find the closest depths
get_closest_depths <- function(df, target_depths) {
  closest_values <- target_depths %>% 
    lapply(function(target) {
      df %>%
        mutate(Distance = abs(Depth - target)) %>%  # Calculating absolute difference
        slice_min(Distance, n = 1) %>%  # Selecting the row with the smallest difference
        mutate(Target_Depth = target)  # Adding a column for target depth
    }) %>%
    bind_rows()  
  return(closest_values)
}

# Applying the function grouped by station
closest_depths_df <- ctd_profiles %>%
  group_by(Location_Station_Number) %>%
  group_modify(~ get_closest_depths(.x, target_depths)) %>%
  ungroup() %>%
  select(Location_Station_Number, Target_Depth, Depth, everything()) %>%
  rename(Measured_depth = Depth)

# After examining the distance between target depth and measured depth, it is evident that both PE477_5 and PE486_3 are not 30 m deep. 
# Also in cases there are more than one values for the same measured_depths, I will take the first one
closest_depths_df <- closest_depths_df %>%
  dplyr::filter(Distance < 2) %>%
  group_by(Location_Station_Number, Target_Depth) %>%
  slice_head(n = 1) %>%  # Select the first row in each group
  ungroup()

write.csv(closest_depths_df, "./results/ctd_profiles/ctd_variables_sampling_depths.csv", row.names = F)
