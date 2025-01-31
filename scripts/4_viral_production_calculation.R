#AIM: To run viralprod to extract lytic and lysogenic  viral production rates 

# 0.0 Setting up ####

source("./scripts/0_source.R")

# 1.0 Importing counts #####

counts_per_mL<- read.csv("results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv")
str(counts_per_mL)

# 2.0 Visual inspection to remove outliers ####
#Plotting VDC against time per assay to check for any outliers ####
# Plotting
vdc_plot <- ggplot(counts_per_mL, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station, scales = "free_y") +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(title = "Loss corrected Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "replicate")
vdc_plot
ggsave(vdc_plot, filename = "./figures/loss_corrected_vdc_plot.svg", width = 20, height = 5)

# SAVING PLOTS PER STATION ####
{
  # Ensure the output folder exists
  output_folder <- "./figures/station_vdc_BEFORE_OUTLIERS_REMOVED_plots"
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Loop through each unique Location_Station
  unique_stations <- unique(counts_per_mL$Location_Station)
  
  for (station in unique_stations) {
    # Filter data for the current station
    station_data <- counts_per_mL %>% dplyr::filter(Location_Station == station)
    
    # Create the plot
    vdc_plot <- ggplot(station_data, aes(x = Timepoint, y = c_Viruses / 1e+6)) +
      geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each replicate
      geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Replicate
      facet_grid(. ~ Sample_Type, scales = "free_y") +  # Facet Sample_Type on left and right
      scale_color_npg() +
      theme_bw() +
      labs(title = paste("Loss corrected Viral Abundance (c_Viruses) - Station", station),
           x = "Timepoint",
           y = "Viral Count (c_Viruses) (millions)",
           color = "Replicate")
    
    # Save the plot
    filename <- paste0(output_folder, "/loss_corrected_vdc_plot_", station, ".svg")
    ggsave(vdc_plot, filename = filename, width = 12, height = 5)
  }
}


# Function to make an exclusion dataframe
add_to_exclusion <- function(existing_df = NULL, location, station, sample_type, timepoint, replicate) {
  
  # Check if the exclusion data frame exists, if not, create it with the proper column names
  if (is.null(existing_df)) {
    existing_df <- data.frame(
      Location = character(),
      station_Number = numeric(),
      Sample_Type = character(),
      Timepoint = numeric(),
      Replicate = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Append the new exclusion row
  new_row <- data.frame(
    Location = location,
    Station_Number = station,
    Sample_Type = sample_type,
    Timepoint = timepoint,
    Replicate = replicate,
    stringsAsFactors = FALSE
  )
  
  # Combine the existing data frame with the new row
  existing_df <- rbind(existing_df, new_row)
  
  return(existing_df)
}


# Initialize the exclusion data frame
# to_exclude <- NULL
# to_exclude <- add_to_exclusion(to_exclude, "PE477", 4, "VPC", 6, 3) # Location, Station_Number, Sampe_Type, Timepoint, Replicate# recursively add on to the dataframe

# write.csv(to_exclude, file = "./results/viral_production_analyses/excluded_samples.csv", row.names =F)

# # To remove in case
# to_exclude <- to_exclude %>%
#   filter(!(Location == "PE477" & 
#              Station_Number == 5 & 
#              Sample_Type == "VPC" & 
#              Timepoint == 24 & 
#              Replicate == 3))

# I have created an exclusion datframe already, and would directly use it to exclude outliers
exclusion_df <- read.csv("./results/viral_production_analyses/excluded_samples_CB.csv")

for (i in 1:nrow(exclusion_df)) {
  counts_per_mL <- counts_per_mL %>%
    mutate(c_Viruses = ifelse(Location == exclusion_df$Location[i] &
                                Station_Number == exclusion_df$Station_Number[i] &
                                Timepoint == exclusion_df$Timepoint[i] &
                                Replicate == exclusion_df$Replicate[i] &
                                Sample_Type == exclusion_df$Sample_Type[i], 
                              NA, c_Viruses))
}


filtered_data <- counts_per_mL %>%
  dplyr::filter(Sample_Type %in% c("VP", "VPC"))  %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

write.csv(filtered_data, "results/viral_loss_corrected_counts/FILTERED_PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv", row.names = F)

# Plotting
vdc_filtered_plot <- ggplot(filtered_data, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
    x = "Timepoint",
    y = "Viral Count (c_Viruses)",
    color = "Replicate")
vdc_filtered_plot
ggsave(vdc_filtered_plot, filename = "./figures/loss_corrected_vdc_filtered_plot.svg", width = 20, height = 5)

vdc_filtered_lm_plot <- ggplot(filtered_data, aes(x = Timepoint, y = c_Viruses / 1e+6)) +
  #geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Replicate
  geom_smooth(aes(group = as.factor(Replicate), color = as.factor(Replicate)), 
              method = "lm", se = FALSE, linetype = "dashed") +  # Add regression lines per replicate
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
    x = "Timepoint",
    y = "Viral Count (c_Viruses)",
    color = "Replicate")
vdc_filtered_lm_plot
ggsave(vdc_filtered_lm_plot, filename = "./figures/loss_corrected_vdc_filtered_lm_plot.svg", width = 20, height = 5)

# ggplot(filtered_data, aes(x = Timepoint, y = c_Viruses / 1e+6)) +
#   #geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
#   geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Replicate
#   geom_smooth(aes(group = as.factor(Replicate), color = as.factor(Replicate)), 
#               method = "lm", se = FALSE, linetype = "dashed") +  # Add regression lines per replicate
#   facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
#   scale_color_npg() +
#   theme_bw(base_size = 15) +
#   labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
#     x = "Timepoint",
#     y = "Viral Count (c_Viruses)",
#     color = "Replicate")


# 3.0 Plotting averaged counts for VP and Diff treatments ####
# using internal viralprod functions available in 0_source.R 


filtered_viruses <- vp_average_replicate_dataframe(filtered_data, add_timepoints = FALSE) %>%
  dplyr::filter(Microbe == "Viruses")

# Creating VP treatment plots
vp_plot <- ggplot(filtered_viruses %>%
         dplyr::filter(Sample_Type == "VP"), 
       aes(x = Timepoint, y = Mean/1e+6, color = Sample_Type)) +
  geom_line(size = 1, color = "#d43028") +  # Line for mean
  geom_point(size = 2, color = "#d43028") +  # Points for mean
  geom_errorbar(aes(ymin = (Mean - SE)/1e+6, ymax = (Mean + SE)/1e+6), width = 0.2, color = "#d43028") +  # Error bars
  facet_wrap(~tag, scales = "free") +  # Facet by tag
  theme_minimal() +
  labs(
    title = "VP treatment for lytic viral production - Mean and SE of loss corrected (T0/T24) viral counts",
    x = "Timepoint",
    y = "Mean viral count (VLPs in millions per mL)",
    color = "Sample Type"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
vp_plot 
ggsave(vp_plot, filename = "./figures/VP_treament_loss_corrected_VDC_filtered_plot.svg", width = 20, height = 10)




# Creating VP treatment plots
diff_plot <- ggplot(filtered_viruses %>%
         dplyr::filter(Sample_Type == "Diff"), 
       aes(x = Timepoint, y = Mean/1e+6, color = Sample_Type)) +
  geom_line(size = 1, color = "#4d778b") +  # Line for mean
  geom_point(size = 2, color = "#4d778b") +  # Points for mean
  geom_errorbar(aes(ymin = (Mean - SE)/1e+6, ymax = (Mean + SE)/1e+6), width = 0.2, color = "#4d778b") +  # Error bars
  facet_wrap(~tag, scales = "free") +  # Facet by tag
  theme_minimal() +
  labs(
    title = "VPC-VP treatment for lysogenic viral production - Mean and SE of loss corrected (T0/T24) viral counts",
    x = "Timepoint",
    y = "Mean viral count difference (in millions per mL)",
    color = "Sample Type"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

diff_plot 
ggsave(diff_plot, filename = "./figures/Diff_treament_loss_corrected_VDC_filtered_plot.svg", width = 20, height = 10)



# 4.0 Running viralprod ####

#Import abundance data
url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
str(abundance)

#Extracting first entries per Location/Station_Number combinations, as some stations had multiple depths and this information is missing from the fiel on GitHub.
#The first one is the depth for VP assays.
abundance<- abundance %>%
  dplyr::filter(Depth == 7,
                Location %in%c("PE477", "PE486"))

#Some checks
vp_check_populations(filtered_data)

class(abundance)
abundance <- vp_class_ori_abu(abundance)
class(abundance)

class(filtered_data) #failed
filtered_data <- vp_class_count_data(filtered_data) #passed
class(filtered_data)

#Running viralprod
 
vp_end_to_end(data = filtered_data ,
              original_abundances = abundance,
              methods = c(2,9,10),
              write_output = T,
              output_dir = 'results/viral_production_analyses/PE_Cruises_0.22_corrected_viral_production_CB_outliers')

vp <- read.csv("results/viral_production_analyses/PE_Cruises_0.22_corrected_viral_production_CB_outliers/vp_results_BP.csv")

vp <- vp %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Sample_Type != "VPC")

vp$Location_Station <- paste(vp$Location, vp$Station_Number, sep = "_")


ggplot(vp, aes(x = Location_Station, y = VP)) +
 # geom_point() +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample_Type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "VP vs Location_Station",
    x = "Location and Station",
    y = "Viral production rate (VLPs per mL per hour)"
  )

# 4.1 Runnning viralprod on T0_T12 corrected couns #####

counts_per_mL_0_12<- read.csv("results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_12_per_mL.csv")
str(counts_per_mL_0_12)
ggplot(counts_per_mL_0_12, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station, scales = "free_y") +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(title = "Loss corrected Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "replicate")
vdc_plot

exclusion_df <- read.csv("./results/viral_production_analyses/excluded_samples.csv")

for (i in 1:nrow(exclusion_df)) {
  counts_per_mL_0_12 <- counts_per_mL_0_12 %>%
    mutate(c_Viruses = ifelse(Location == exclusion_df$Location[i] &
                                Station_Number == exclusion_df$Station_Number[i] &
                                Timepoint == exclusion_df$Timepoint[i] &
                                Replicate == exclusion_df$Replicate[i] &
                                Sample_Type == exclusion_df$Sample_Type[i], 
                              NA, c_Viruses))
}


filtered_data_0_12 <- counts_per_mL_0_12 %>%
  dplyr::filter(Sample_Type %in% c("VP", "VPC"))  %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

# Plotting
ggplot(filtered_data_0_12, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
    x = "Timepoint",
    y = "Viral Count (c_Viruses)",
    color = "Replicate")


vp_check_populations(filtered_data_0_12)

class(abundance)
abundance <- vp_class_ori_abu(abundance)
class(abundance)

class(filtered_data_0_12) #failed
filtered_data_0_12 <- vp_class_count_data(filtered_data_0_12) #passed
class(filtered_data_0_12)



vp_end_to_end(data = filtered_data_0_12 ,
              original_abundances = abundance,
              methods = c(2,9,10),
              write_output = T,
              output_dir = 'results/viral_production_analyses/PE_Cruises_0.22_0_12corrected_viral_production')

vp_0_12 <- read.csv("results/viral_production_analyses/PE_Cruises_0.22_0_12corrected_viral_production/vp_results_BP.csv")

vp_0_12 <- vp_0_12 %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Sample_Type != "VPC")

vp_0_12$Location_Station <- paste(vp_0_12$Location, vp_0_12$Station_Number, sep = "_")


ggplot(vp_0_12, aes(x = Location_Station, y = VP)) +
  # geom_point() +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample_Type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "VP vs Location_Station",
    x = "Location and Station",
    y = "VP"
  )


# hAVE to correct it for bacterail eficiency t sho it coria ####
vp_0_12_wide <- vp_0_12 %>%
  select(c(Location_Station, Sample_Type, VP)) %>%
  pivot_wider(names_from = Sample_Type,
              values_from = VP) 
bac_eff <- vp %>%
  select(c(Location, Station_Number, bac_efficiency))%>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

vp_0_12_wide <- vp_0_12_wide %>%
  left_join(bac_eff, by = "Location_Station") %>%
  mutate(corr_VP = VP*100/bac_efficiency,
         corr_Diff = Diff*100/bac_efficiency)

### combining vp  and vp_0_12


# Add a grouping column to distinguish datasets
vp <- vp %>% mutate(Group = "vp_0_24")
vp_0_12 <- vp_0_12 %>% mutate(Group = "vp_0_12")

# Filter `vp` to retain only the `Location_Station` values from `vp_0_12`
vp_filtered <- vp %>% dplyr::filter(Location_Station %in% vp_0_12$Location_Station)

# Combine the two datasets
combined_data <- bind_rows(vp_filtered, vp_0_12)

# Create the grouped bar plot
ggplot(combined_data, aes(x = Location_Station, y = VP, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bar plot
  facet_wrap(~ Sample_Type) +  # Facet by Sample_Type
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "VP vs Location_Station (Grouped by vp_0_24 and vp_0_12)",
    x = "Location and Station",
    y = "VP",
    fill = "Group"
  )
vp_data_loss <- combined_data %>% dplyr::filter(Sample_Type == "VP")
kruskal.test(VP ~ Group, data = vp_data_loss)

diff_data_loss <- combined_data %>% dplyr::filter(Sample_Type == "Diff")
kruskal.test(VP ~ Group, data = diff_data_loss)


# Coapring lytic viral production andlysogenic viralproduction betwen T0_T12 an T0_T24 #######
vp_all<- read.csv("results/viral_production_analyses/PE_Cruises_0.22_corrected_viral_production_CB_outliers/vp_results_ALL.csv") %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE", 
                Sample_Type != "VPC",
                Time_Range %in% c("T0_T12", "T0_T24"))
vp_all$Location_Station <- paste(vp_all$Location, vp_all$Station_Number, sep = "_")

ggplot(vp_all, aes(x = Location_Station, y = VP, fill = Time_Range)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bar plot
  facet_wrap(~ Sample_Type) +  # Facet by Sample_Type
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text
    strip.text = element_text(size = 12, face = "bold")  # Style facet labels
  ) +
  labs(
    title = "VP vs Location_Station (T12 vs T24)",
    x = "Location and Station",
    y = "VP",
    fill = "Time Range"
  )


vp_all_wide <- vp_all %>%
  select(c(Location_Station, Sample_Type, Time_Range, VP)) %>%
  pivot_wider(names_from = Sample_Type,
              values_from = VP) 
bac_eff <- readxl::read_excel("metadata/PE477_PE486_VP_LM_HMS.xlsx", sheet = 'final_output') %>%
  select(c(Location, Station_Number, bac_efficiency))%>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

vp_all_wide <- vp_all_wide %>%
  left_join(bac_eff, by = "Location_Station") %>%
  mutate(corr_Lytic = VP*100/bac_efficiency,
         corr_Lysogenic = Diff*100/bac_efficiency)

ggplot(vp_all_wide, aes(x = Location_Station, y = corr_Lytic, fill = Time_Range)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bar plot
 # facet_wrap(~ Sample_Type) +  # Facet by Sample_Type
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text
    strip.text = element_text(size = 12, face = "bold")  # Style facet labels
  ) +
  scale_fill_futurama()+
  labs(
    title = "Lytic viral production rate (T12 vs T24)",
    x = "Location and Station",
    y = "Viral production (VLPs per mL per hour)",
    fill = "Time Range"
  )

ggplot(vp_all_wide, aes(x = Location_Station, y = corr_Lysogenic, fill = Time_Range)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bar plot
  # facet_wrap(~ Sample_Type) +  # Facet by Sample_Type
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text
    strip.text = element_text(size = 12, face = "bold")  # Style facet labels
  ) +
  scale_fill_futurama()+
  labs(
    title = "Lysogenic viral production rate (T12 vs T24)",
    x = "Location and Station",
    y = "Viral production (VLPs per mL per hour)",
    fill = "Time Range"
  )

vp_data <- vp_all %>% dplyr::filter(Sample_Type == "VP")
kruskal.test(VP ~ Time_Range, data = vp_data)

diff_data <- vp_all %>% dplyr::filter(Sample_Type == "Diff")
kruskal.test(VP ~ Time_Range, data = diff_data)




# Bacterial efficiency corrceted 

bac_eff_df <- read.csv("./results/PE477_PE486_3depths_combined_7m.csv") %>%
  dplyr::select(c("Location", "Station_Number", "bac_efficiency"))
vp <- vp %>%
  left_join(bac_eff_df, by = c("Location", "Station_Number")) %>%
  mutate(eff_corrected_VP = (VP *100)/ bac_efficiency)

ggplot(vp, aes(x = Location_Station, y = eff_corrected_VP)) +
  # geom_point() +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample_Type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "bacterial loss corrected VP vs Location_Station",
    x = "Location and Station",
    y = "VP"
  )



# Viraprod wihtout removing outliers ####

counts_per_mL<- read.csv("results/viral_loss_corrected_counts/PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv")
str(counts_per_mL)



class(counts_per_mL) #failed
counts_per_mL <- vp_class_count_data(counts_per_mL) #passed
class(counts_per_mL)

#Running viralprod

vp_end_to_end(data = counts_per_mL ,
              original_abundances = abundance,
              methods = c(2,9,10),
              write_output = T,
              output_dir = 'results/viral_production_analyses/PE_Cruises_0.22_corrected_with_outliers_viral_production')

vp_with_outliers <- read.csv("results/viral_production_analyses/PE_Cruises_0.22_corrected_with_outliers_viral_production/vp_results_BP.csv")

vp_with_outliers <- vp_with_outliers %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Sample_Type != "VPC")

vp_with_outliers$Location_Station <- paste(vp_with_outliers$Location, vp_with_outliers$Station_Number, sep = "_")


ggplot(vp_with_outliers, aes(x = Location_Station, y = VP)) +
  # geom_point() +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample_Type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "VP vs Location_Station",
    x = "Location and Station",
    y = "VP"
  )


# Combined plot of outliers removed with ones not removed

vp <- vp %>% mutate(Data_Source = "VP_without_outliers")
vp_with_outliers <- vp_with_outliers %>% mutate(Data_Source = "VP_with_Outliers")

# Combining both dataframes
vp_combined <- bind_rows(vp, vp_with_outliers)

# Creating the grouped bar plot
ggplot(vp_combined, aes(x = Location_Station, y = VP, fill = Data_Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Sample_Type) +
  labs(title = "Comparison of VP with and without outliers",
       x = "Location Station",
       y = "VP",
       fill = "Data Source") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Total bactera lysed ####

#Assuming a burst size of 50

ggplot(vp, aes(x = Location_Station, y = eff_corrected_VP/50)) +
  # geom_point() +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sample_Type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "bacterial loss rate (lytic/lysogeny) vs Location_Station",
    x = "Location and Station",
    y = "bacterial loss rate in original sample"
  )
