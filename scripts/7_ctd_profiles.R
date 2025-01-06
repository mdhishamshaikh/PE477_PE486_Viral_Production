# AIM: Extract data from CTD profiles

# 0.0 Setting up ####
library(tidyverse)
library(ggsci)

# 1.0 Importing CTD profiles and selecting stations ####
file_path <- "./Data/ctd_profiles"
file_list <- list.files(path = file_path, pattern = "*.csv", full.names = TRUE)

# Function to add file name as a column
read_and_label <- function(file) {
 
  df <- read.csv(file)
  df$file_name <- basename(file)
  return(df)
}


# Changing file names and keeping nly the staions needed.
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
# removing station PE477_7 and PE486_8 as tehy don't have viral production data attached to them.

combined_df <- combined_df %>%
  dplyr::filter(!Location_Station_Number %in% c("PE477_7", "PE486_8")) 

# 2.0 Selecting variables and simplifying names ####

combined_df <- combined_df %>%
  rename(
    Depth = depSM,
    Salinity = sal00,
    Temperature = t090C,
    Pressure = prDM,
    Conductivity = C0,
    Oxygen = sbeox0Mm.L,
    Turbidity = turb,
    Fluorescence = flC
  )

# 3.0 Calculating density and removing surface values values ####

# Gibbs SeaWater 
combined_df$Density <- gsw::gsw_rho(combined_df$Salinity, combined_df$Temperature, combined_df$Pressure)

combined_df <- combined_df %>%
  dplyr::filter(Depth > 3) # Making sure inital first half meter is not considereed
# 4.0 Reshaping and plotting#####

ctd_profiles<- combined_df %>%
  dplyr::select(c("Location_Station_Number", "Depth", "Salinity", "Temperature", "Pressure", "Density", "Conductivity", "Oxygen", "Turbidity", "Fluorescence")) %>%
  separate(Location_Station_Number, c("Location", "Station_Number"), "_", remove = F)

write.csv(ctd_profiles, "./results/ctd_profiles/ctd_profiles.csv", row.names = F)


ctd_long_df <- ctd_profiles %>%
  dplyr::select(-Location, -Station_Number) %>%
  reshape2::melt(id.vars = c("Location_Station_Number", "Depth")) %>%
  mutate(variable = factor(variable, levels = c(
    "Salinity", "Temperature", "Pressure", "Density", 
    "Conductivity", "Oxygen", "Turbidity", "Fluorescence"
  )))

unique(ctd_long_df$variable)

# Plotting
pdf("./results/ctd_profiles/depth_profiles_by_station.pdf", width = 10, height = 7)


for (station in unique(ctd_long_df$Location_Station_Number)) {
  
 
  station_df <- subset(ctd_long_df, Location_Station_Number == station)
 
  p <- ggplot(station_df, aes(x = value, y = Depth)) +
    geom_point(size = 0.1) +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +  # Facet by parameter (variable)
    labs(title = paste("Depth Profiles for", station),
         x = "Parameter Value",
         y = "Depth (m)") +
    scale_y_reverse() +  # Reverse the y-axis for depth profiles
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
          strip.text = element_text(size = 8))              # Adjust facet label size
  
  print(p)
}

dev.off()

# As images

# Create a folder for the images
dir.create("./results/ctd_profiles/per_station", showWarnings = FALSE)

for (station in unique(ctd_long_df$Location_Station_Number)) {
  
  # Subset the data for the current station
  station_df <- subset(ctd_long_df, Location_Station_Number == station)
  
  # Create the plot
  p <- ggplot(station_df, aes(x = value, y = Depth)) +
    geom_point(size = 0.1) +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +  # Facet by parameter (variable)
    labs(title = paste("Depth Profiles for", station),
         x = "Parameter Value",
         y = "Depth (m)") +
    scale_y_reverse() +  # Reverse the y-axis for depth profiles
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
          strip.text = element_text(size = 8))              # Adjust facet label size
  
  # Save the plot as an image
  ggsave(
    filename = paste0("./results/ctd_profiles/per_station/depth_profile_", station, ".png"),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300  # Set a high resolution
  )
}



# Use this to figure ut if clines exist per station and start and end, maybe an avergae too.

# 5.0 Extracing CTTD variables for 7, 15 and 30 m ####

# Target depths
target_depths <- c(7, 15, 30) # Some station might not have 30 m. Use the distance variable to fgure out. Maybe worth adding max depth

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
  rename(Measured_depth = Depth)# Organize columns

# After examining the distance between target depth and measured depth, it is evident that both PE477_5 and PE486_3 are not as deep as 30 m. 
# Also incase there are more than one values for the same measured_depths, I will take the first one
closest_depths_df <- closest_depths_df %>%
  dplyr::filter(Distance < 2)

write.csv(closest_depths_df, "./results/ctd_profiles/ctd_variables_sampling_depths.csv", row.names = F)



# 6.0 CTD clustering ######

library(vegan)
library(ggfortify)
library(psych)

# Normalize CTD variables
ctd_scaled <- ctd_profiles%>%
  select(-Location_Station_Number, -Location, -Station_Number, -Depth, -Pressure) %>%
  scale()

# Performing factor analysis to remove redundancy




# Distance matrix (Euclidean or Bray-Curtis)
dist_matrix <- dist(ctd_scaled, method = "euclidean")

# Hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Assign clusters (e.g., 3 clusters)
clusters <- cutree(hclust_result, k = 3)

# Add cluster information to CTD data
ctd_data <- ctd_profiles %>%
  select(-Pressure) %>%
  mutate(Cluster = as.factor(clusters))



# PCA on CTD data
pca_result <- prcomp(ctd_data %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, -Cluster), scale. = TRUE)

# PCA plot
custom_palette <- pal_d3("category20")(16)  # Generate 16 distinct colors from the D3 palette

autoplot(pca_result, data = ctd_data, shape = "Cluster", color = "Location_Station_Number", size = 3) +
  scale_color_manual(values = custom_palette) +  # Use the custom palette
  labs(title = "PCA of CTD Data", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "right")


custom_palette <- pal_d3("category20")(16)

# Extract PCA loadings
loadings <- as.data.frame(pca_result$rotation)
loadings$Variable <- rownames(loadings)  # Add variable names
loadings <- loadings %>% 
  mutate(PC1 = PC1 * 0.05 ,  # Scale for plotting
         PC2 = PC2 * 0.05)

# PCA plot with loadings
autoplot(pca_result, data = ctd_data, shape = "Cluster", color = "Location_Station_Number", size = 3) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", size = 0.8) +
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = Variable), 
            color = "black", size = 4, vjust = 1.2, hjust = -.2) +
  scale_color_manual(values = custom_palette) +  # Use the custom palette
  labs(title = "PCA of CTD Data with Loadings", x = "PC1", y = "PC2") +
  xlim()
  theme_minimal() +
  theme(legend.position = "right")


# I dentifying split stations

# Summarize clusters by station and depth
cluster_summary <- ctd_data %>%
  group_by(Location_Station_Number, Depth, Cluster) %>%
  summarise(Count = n(), .groups = "drop")

# View summary
print(cluster_summary)
# Count the number of clusters per station and depth
split_stations <- cluster_summary %>%
  group_by(Location_Station_Number, Depth) %>%
  summarise(Num_Clusters = n_distinct(Cluster), .groups = "drop") %>%
  filter(Num_Clusters > 1)  # Keep only those with multiple clusters

# View split stations
print(split_stations)


# 6.2   Clustering only PE477 stations

pe477 <- ctd_profiles %>%
  dplyr::filter(Location == "PE477")

# Factor analysis to remove redundancies

