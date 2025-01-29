# AIM:; To correct FCS counts (ONLY total viruses) for viral loss using 0.22 µm treatments

# 0.0 Setting up ####
source("scripts/0_source.R")

# 1.0 Importing FCS per mL counts for VP, VPC and 0.22 µm treatments ####
#These are TE corrected counts (per mL)
counts_per_mL<- read.csv("PE_Cruises_FCS_with_0.22/results/PE_Cruises_FCS_with_0.22_per_mL.csv")
# Creating a combined variable for Location and Station_Number
counts_per_mL$Location_Station <- paste(counts_per_mL$Location, counts_per_mL$Station_Number, sep = "_")


# 2.0 Visualizing 0.22 µm counts over time ####

# removing outliers
counts_per_mL_outliers_removed <- counts_per_mL %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE477" & Station_Number == 2 & Timepoint == 6 & Replicate == 1,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  )) %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE477" & Station_Number == 6 & Timepoint == 0 & Replicate == 2,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  )) %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE486" & Station_Number == 5 & Timepoint == 9 & Replicate == 3,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  ))



# Create the plot
ggplot(counts_per_mL_outliers_removed %>% dplyr::filter(Sample_Type=="0.22"), aes(x = Timepoint, y = c_Viruses, color = as.factor(Replicate))) +
  # geom_line(aes(group = interaction(Location_Station, Replicate)), size = 1) +
  geom_smooth(method = 'lm', se= F) +
  geom_point(size = 3) +
  facet_wrap(~ Location_Station, scales = "free") +
  labs(
    title = "0.22 µm Viral Counts Over Time",
    x = "Timepoint",
    y = "Viral Counts (c_Viruses)",
    color = "Replicate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )


# 3.0 Calculating slopes to correct VP and VPC counts using 0.22 µm treatment ####
#After having tried a few ways the best way to deal with them is to only take T0 and T24

# I will calculate the slope using these two time points, and add the counts to each time point (except T0)
all_stations_0_24_0.22_slopes<- counts_per_mL_outliers_removed %>%
  dplyr::filter(Sample_Type == "0.22",
                Timepoint %in% c(0,24)) %>%
  group_by(Location_Station) %>%
  summarise(Slope = coef(lm(c_Viruses ~ Timepoint))[2]) %>%
  expand_grid(Timepoint = c(3,6,9,12,24))


time_elapsed<- counts_per_mL_outliers_removed %>%
  dplyr::select(Location_Station, Timepoint) %>%
  distinct() %>%
  arrange(Location_Station, Timepoint) %>%
  group_by(Location_Station) %>%
  mutate(Time_Elapsed = c(0, diff(Timepoint)))


viral_count_correction_per_timepoint_0_24 <- time_elapsed %>%
  left_join(all_stations_0_24_0.22_slopes, by = c("Location_Station", "Timepoint")) %>%
  dplyr::mutate(viral_count_correction = Time_Elapsed * abs(Slope))


loss_corrected_viral_counts_0_24 <- counts_per_mL_outliers_removed %>%
  dplyr::filter(Sample_Type != "0.22") %>%
  left_join(
    viral_count_correction_per_timepoint_0_24 %>%
      dplyr::mutate(viral_count_correction = if_else(Timepoint == 0, 0, viral_count_correction)),
    by = c("Location_Station", "Timepoint")
  )  %>%
  dplyr::mutate(loss_corrected_c_Viruses = c_Viruses + viral_count_correction) %>%
  dplyr::mutate(before_correction_c_Viruses = c_Viruses) %>%
  dplyr::mutate(c_Viruses = loss_corrected_c_Viruses) %>%
  dplyr::select(-c("c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3", "VBR", "HNAperLNA"))
  

write.csv(loss_corrected_viral_counts_0_24, "results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv")


# Visualizing all stations corrected with T0 and T24 0.22 µm slopes
counts_per_mL_0_24<- read.csv("results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv")
vdc_plot_0_24 <- ggplot(counts_per_mL_0_24  %>%
  dplyr::filter(Sample_Type != "0.22") , aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(title = "Loss corrected (T0 & T24) Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "replicate")
vdc_plot_0_24

vdc_plot_0_24_high24_stations <- ggplot(counts_per_mL_0_24  %>%
                          dplyr::filter(Sample_Type != "0.22",
                                        Location_Station %in% c("PE477_2", "PE477_4", "PE486_2", "PE486_4", "PE486_6", "PE486_7")) , 
                          aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  ylim(c(0,13)) +
  theme_bw(base_size = 15) +
  labs(title = "Stations with high T24 - Loss corrected (T0 & T24) Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "replicate")
vdc_plot_0_24_high24_stations


# Correcting the Stations with high T24  with only T0 to T12 slopes

high_24_stations_0_12_0.22_slopes<- counts_per_mL_outliers_removed %>%
  dplyr::filter(Sample_Type == "0.22",
                Timepoint %in% c(0, 3, 6, 9, 12),
                Location_Station %in% c("PE477_2", "PE477_4", "PE486_2", "PE486_4", "PE486_6", "PE486_7")) %>%
  group_by(Location_Station) %>%
  summarise(Slope = coef(lm(c_Viruses ~ Timepoint))[2]) %>%
  expand_grid(Timepoint = c(3,6,9,12,24))

viral_count_correction_per_timepoint_0_12 <- time_elapsed %>%
  left_join(high_24_stations_0_12_0.22_slopes, by = c("Location_Station", "Timepoint")) %>%
  dplyr::mutate(viral_count_correction = Time_Elapsed * abs(Slope))


loss_corrected_viral_counts_0_12 <- counts_per_mL_outliers_removed %>%
  dplyr::filter(Sample_Type != "0.22",
                Location_Station %in% c("PE477_2", "PE477_4", "PE486_2", "PE486_4", "PE486_6", "PE486_7")) %>%
  left_join(
    viral_count_correction_per_timepoint_0_12 %>%
      dplyr::mutate(viral_count_correction = if_else(Timepoint == 0, 0, viral_count_correction)),
    by = c("Location_Station", "Timepoint")
  )  %>%
  dplyr::mutate(loss_corrected_c_Viruses = c_Viruses + viral_count_correction)%>%
  dplyr::mutate(before_correction_c_Viruses = c_Viruses) %>%
  dplyr::mutate(c_Viruses = loss_corrected_c_Viruses) %>%
  dplyr::select(-c("c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3", "VBR", "HNAperLNA"))

write.csv(loss_corrected_viral_counts_0_12, "results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_12_per_mL.csv")

counts_per_mL_0_12<- read.csv("results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_12_per_mL.csv")
vdc_plot_0_12_high24_stations <- ggplot(counts_per_mL_0_12, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  ylim(c(0,13)) +
  theme_bw(base_size = 15) +
  labs(title = "Stations with high T24 - Loss corrected (T0-T12) Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "replicate")
vdc_plot_0_12_high24_stations

