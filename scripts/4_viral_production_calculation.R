#AIM: To run viralprod to extract lytic and lysogenic  viral production rates 

# 0.0 Setting up ####

source("./scripts/0_source.R")

# 1.0 Importing counts #####

counts_per_mL<- read.csv("./PE_Cruises_FCS/results/PE_Cruises_FCS_per_mL.csv")
str(counts_per_mL)
counts_per_mL<- counts_per_mL %>%
  dplyr::filter(counts_per_mL$Sample_Type != '0.22') # 0.22 will be used for viral decay

#There's some additional files that I would like to exclude from this. These are replicates with low events/sec
#To do so I will use the selected data from the metadata file. This should have already been done before processing. 
#Importing the .xlsx file

selected_files<- read_excel("Metadata/PE_Cruises_VP_Metadata.xlsx", sheet = 'Selected_Metadata') #contains selected filenames

counts_per_mL<- inner_join(counts_per_mL, selected_files, by = "Sample_Name")
#write.csv(counts_per_mL, file = 'Linear_Regression_Method_Supervised/PE477_PE486_filtered.csv', row.names = F)

#Plotting VSC against time per assay to check for any outliers ####

filtered_data <- counts_per_mL %>%
  dplyr::filter(Sample_Type %in% c("VP", "VPC"))  %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

# Plotting
vdc_plot <- ggplot(filtered_data, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station) +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
       x = "Timepoint",
       y = "Viral Count (c_Viruses)",
       color = "Sample Type")
vdc_plot
ggsave(vdc_plot, filename = "./figures/vdc_plot.svg", width = 20, height = 5)

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
exclusion_df <- read.csv("./results/viral_production_analyses/excluded_samples.csv")

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
ggsave(vdc_filtered_plot, filename = "./figures/vdc_filtered_plot.svg", width = 20, height = 5)

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
ggsave(vdc_filtered_lm_plot, filename = "./figures/vdc_filtered_lm_plot.svg", width = 20, height = 5)

##########################
library(dplyr)
library(ggplot2)

# Step 1: Group by Location, Station_Number, and Timepoint, then average replicates for each treatment
averaged_data <- filtered_data %>%
  group_by(Location, Station_Number, Timepoint, Sample_Type) %>%
  summarise(mean_c_Viruses = mean(c_Viruses, na.rm = TRUE),
            sem_c_Viruses = sd(c_Viruses, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

# Step 2: Separate VP and VPC treatments and calculate the difference
wide_data <- averaged_data %>%
  pivot_wider(names_from = Sample_Type, 
              values_from = c(mean_c_Viruses, sem_c_Viruses), 
              names_sep = "_")

# Step 3: Calculate the difference between VPC and VP treatments and their standard error
wide_data <- wide_data %>%
  mutate(
    difference_c_Viruses = mean_c_Viruses_VPC - mean_c_Viruses_VP,  # Difference between VPC and VP
    difference_sem = sqrt((sem_c_Viruses_VPC)^2 + (sem_c_Viruses_VP)^2)  # Simplified error calculation
  )

# Step 4: Plot the difference curve with standard error
ggplot(wide_data, aes(x = Timepoint, y = difference_c_Viruses / 1e+6)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = (difference_c_Viruses - difference_sem) / 1e+6,
                    ymax = (difference_c_Viruses + difference_sem) / 1e+6), 
                width = 0.2) +
  facet_wrap(~ paste(Location, Station_Number, sep = "_"), scales = "free") +
  theme_bw(base_size = 15) +
  labs(title = "Difference Curve: VPC - VP",
       x = "Timepoint",
       y = "Difference in Viral Count (c_Viruses in millions)")
###########

# 
# 
# wide_data
# 
# wide_data_PE477_3 <- wide_data %>%
#   dplyr::filter(Location == "PE477",
#                 Station_Number == 3) %>%
#   arrange(Timepoint)  # Ensure timepoints are in order
# 
# wide_data_PE477_1 <- wide_data %>%
#   dplyr::filter(Location == "PE477",
#                 Station_Number == 1) %>%
#   arrange(Timepoint)  # Ensure timepoints are in order
# 
# 
# DF2 <-wide_data_PE477_3
# DF2 <-wide_data_PE477_1
# # Sample DataFrame to test the function
# index_peaks <- vp_determine_peaks(c(+10e+10, DF2$difference_c_Viruses , -10e+10))
# index_valleys <- vp_determine_valleys(c(+10e+10, DF2$difference_c_Viruses , -10e+10))
# index_peaks
# index_valleys
# # Peak at 5 and valley at 3
# 
# index_peaks <- vp_determine_peaks_with_se(c(+10e+10, DF2$difference_c_Viruses, -10e+10),
#                                           c(0, DF2$difference_sem , 0))
# index_valleys <- vp_determine_valleys_with_se(c(+10e+10, DF2$difference_c_Viruses, -10e+10),
#                                               c(0, DF2$difference_sem , 0))
# 
# index_peaks
# index_valleys
# 
# 
# 
# 
# 
# if (length(index_peaks) == 0){
#   viral_production <- 0
#   abs_vp <- 0
#   se <- 0
# }else {
#   total_vp <- 0
#   total_abs_vp <- 0
#   total_se <- 0
#   
#   for (index in 1:length(index_peaks)){
#     viral_production_index <- (DF2$difference_c_Viruses[index_peaks[index]] - DF2$difference_c_Viruses[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
#     total_vp <- total_vp + viral_production_index
#     
#     abs_vp_index <- DF2$difference_c_Viruses[index_peaks[index]] - DF2$difference_c_Viruses[index_valleys[index]]
#     total_abs_vp <- total_abs_vp + abs_vp_index
#     
#     se_index <- (DF2$difference_sem[index_peaks[index]] + DF2$difference_sem[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
#     total_se <- total_se + se_index
#   }
#   viral_production <- total_vp / length(index_peaks)
#   abs_vp <- total_abs_vp
#   se <- total_se / length(index_peaks)
# }
# 
# 
# result <- c(combi_tag, time, virus, sample, viral_production, abs_vp, se)
# result_list[[length(result_list) + 1]] <- result
# }
# 
# 
# 
# # Display the sample dataframe
# print(AVG_dataframe)
# 
# # Check for absence of peaks and valleys
# if (length(index_peaks) == 0 || length(index_valleys) == 0) {
#   viral_production <- 0
#   abs_vp <- 0
#   se <- 0
# } else {
#   total_vp <- 0
#   total_abs_vp <- 0
#   total_se <- 0
#   
#   # Calculate viral production, absolute VP, and standard error (SE) based on peak-valley pairs
#   for (index in 1:min(length(index_peaks), length(index_valleys))) {
#     # Calculate production index using `difference_c_Viruses`
#     viral_production_index <- (DF2$difference_c_Viruses[index_peaks[index]] - DF2$difference_c_Viruses[index_valleys[index]]) / 
#       (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
#     total_vp <- total_vp + viral_production_index
#     
#     # Calculate absolute VP difference
#     abs_vp_index <- DF2$difference_c_Viruses[index_peaks[index]] - DF2$difference_c_Viruses[index_valleys[index]]
#     total_abs_vp <- total_abs_vp + abs_vp_index
#     
#     # Calculate SE using `difference_sem`
#     se_index <- (DF2$difference_sem[index_peaks[index]] + DF2$difference_sem[index_valleys[index]]) / 
#       (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
#     total_se <- total_se + se_index
#   }
#   
#   # Average the total viral production and SE over the number of peak-valley pairs
#   viral_production <- total_vp / length(index_peaks)
#   abs_vp <- total_abs_vp
#   se <- total_se / length(index_peaks)
# }

##########


#Import abundance data
url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
str(abundance)

#Extracting first entries per Location/Station_Number combinations, as some stations had multiple depths and this information is missing from the fiel on GitHub.
#The first one is the depth for VP assays.
abundance<- abundance %>%
  dplyr::filter(Depth %in% c(1, 7) )

#Some checks
vp_check_populations(counts_per_mL)

class(abundance)
abundance <- vp_class_ori_abu(abundance)
class(abundance)

class(counts_per_mL) #failed
counts_per_mL <- vp_class_count_data(counts_per_mL) #passed
class(counts_per_mL)

#Running viralprod
 
vp_end_to_end(data = counts_per_mL ,
              original_abundances = abundance,
              methods = c(2,9,10),
              write_output = T,
              output_dir = 'results/viral_production_analyses/PE_Cruises_viral_production')


#SOME ADDITIONAL INSPECTION ####
# #### Visualize ####
# 
# cruise_vp<- read.csv("PE_Cruises/results/PE_Cruises_viral_production/vp_results_BP.csv")
# cruise_vp<- cruise_vp %>%
#   dplyr::filter(VP_Method == 'VPCL_AR_DIFF_SE',
#                 Sample_Type != 'VPC',
#                 Population == 'c_Viruses') %>%
#   mutate(Stations = paste0(Location, "_",as.character(Station_Number)))
# str(cruise_vp)
# 
# cruise_vp
# 
# ggplot(data = cruise_vp,
#        aes(x = Stations,
#            y = abs_VP,
#            color = Sample_Type,
#            fill = Time_Range))+
#   geom_point()
# 
# 
# abundance7<- abundance %>%
#   dplyr::filter(Depth %in% c(1,7),
#                 Location != "NJ2020") %>%
#   mutate(Depth = as.integer(Depth)) 
# str(abundance7)
# 
# vp<- vp_results_output_BP_df %>%
#   dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
#                 Sample_Type != 'VPC',
#                 Population == 'c_Viruses') %>%
#   mutate(Sample_Type = base::ifelse(Sample_Type == 'VP', 'Lytic', 'Lysogenic')) %>%
#   select(-c(VP_Method, VP_R_Squared, VP_SE)) %>%
#   pivot_wider(names_from = Sample_Type,
#               values_from = c(VP, abs_VP))
# vp$Depth<- 7
# str(vp)
# #need to wrangle it tolytic andlysogenic prod columns
# 
# vp_abundance <- full_join(vp, abundance7, by = c("Location", "Station_Number", "Depth"))
# head(vp_abundance)
# 
# #REMOVE ROWS THAT DON'T HAVE VP assay perfomred.
# 
# vp_abundance<- vp_abundance %>% 
#   dplyr::filter(!(is.na(VP_Lytic)))
# 
# #This is the combined datframe. Save it as csv and we'll use thisf or visualization.
# 
# write.csv(vp_abundance, paste0(getwd(), 
#                                "/PE_Cruises/results/PE_Cruises_viral_production/vp_abundance_nutrients_tsV.csv"),
#           row.names = F)
# 
# 
# # Combining this with supervised data
# lm_supervised <- read_excel("Linear_Regression_Method_Supervised/PE477_PE486_VP_LM.xlsx", 
#                             sheet = "final_output") %>%
#   rename(LM_Lytic = corrected_Lytic, LM_Lysogenic = corrected_Lysogenic) %>%
#   mutate(station_code = paste(Location, Station_Number, Depth, sep = "_")) %>%
#   select(c(station_code, LM_Lytic, LM_Lysogenic))
# 
# vp_abundance_ts <- vp_abundance %>% 
#   mutate(station_code = paste (Location,Station_Number,Depth, sep = "_"))
# 
# vp_pe_df <- merge(vp_abundance_ts, lm_supervised, by = "station_code" )
# 
# write.csv(vp_pe_df, paste0(getwd(), 
#                                "/PE_Cruises/results/PE_Cruises_viral_production/vp_abundance_nutrients_ts_LM_sup.txt"),
#           row.names = F)
# 
# lm_vs_vpcl <- vp_pe_df %>%
#   select(station_code, VP_Lytic, VP_Lysogenic, LM_Lytic, LM_Lysogenic) %>%
#   pivot_longer(cols = -station_code,
#                names_to = c("analysis_method", "VP_type"), 
#                values_to = 'VP',
#                names_sep = '_') %>%
#   mutate(analysis_method = recode(analysis_method,
#                              'VP' = 'VIPCAL',
#                              'LM' = 'LM-S')) %>%
#   mutate(analysis_method = factor(analysis_method, levels = c("LM-S", "VIPCAL")),
#          VP_type = factor(VP_type, levels = c("Lytic", "Lysogenic"))) %>%
#   mutate(VP = ifelse(VP < 0, 0, VP))
# 
# 
# # Visualize
# 
# ggplot(lm_vs_vpcl ,  
#        aes(x = station_code, y = VP, fill = analysis_method
#            )) +
#   geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
#   # geom_errorbar(aes(ymin = VP - VP_SE, ymax = VP + VP_SE),
#   #               position = position_dodge(width = 0.7), width = 0.25) +
#   labs(title = "PE477 & PE486 Viral Production",
#        x = "Station Code",
#        y = "Viral Production (VP)",
#        fill = "Analytical Method") +
#   scale_fill_manual(values = c("VIPCAL" = "#d43028", "LM-S" = "#4d778b")) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
#   facet_wrap(~ VP_type, ncol = 2, scales = "fixed") +  # Use 'facet_wrap' with free y-axis scaling
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(angle = 90))
# 
# 
# ### Percentage bacteria lysed and percent lysogeny
# # Assuming burst sizes - 10, 25, and 40
# 
# 
# burst_sizes <- c(10, 25, 40)
# 
# calc_lysed_lysogeny <- function(df, burst_sizes) {
#   results <- df %>%
#     select(station_code, Total_Bacteria, starts_with("VP_"), starts_with("LM_")) %>%
#     pivot_longer(cols = starts_with("VP_") | starts_with("LM_"),
#                  names_to = c("analysis_method", "VP_type"),
#                  values_to = "VP",
#                  names_pattern = "(VP|LM)_(Lytic|Lysogenic)") %>%
#     filter(!is.na(Total_Bacteria)) %>%  # Filter out rows with NA in Total_Bacteria
#     rowwise()
#   
#   final_results <- list()
#   
#   for (burst_size in burst_sizes) {
#     temp_results <- results %>%
#       mutate(burst_size = burst_size,
#              bacteria_Lys = VP / burst_size,
#              percent_bacteria_lys = 100 * (VP / burst_size) / Total_Bacteria)
#     
#     final_results <- bind_rows(final_results, temp_results)
#   }
#   
#   final_results
# }
# 
# # Apply the function to the data frame
# vp_pe_df_combined <- calc_lysed_lysogeny(vp_pe_df, burst_sizes)
# 
# # Display the combined dataframe
# print(vp_pe_df_combined)
# 
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# # Ensure the data has a 'bacteria_lysogeny' column calculated from lysogenic types
# vp_pe_df_combined <- vp_pe_df_combined %>%
#   mutate(bacteria_lysogeny = if_else(VP_type == "Lysogenic", bacteria_Lys, NA_real_))
# 
# # Prepare data for plotting
# plot_data <- vp_pe_df_combined %>%
#   filter(VP_type %in% c("Lytic", "Lysogenic")) %>%
#   mutate(station_code = factor(station_code, levels = unique(station_code)))
# 
# # Plotting
# ggplot(plot_data, aes(x = station_code, y = bacteria_Lys, fill = analysis_method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_grid(burst_size ~ VP_type, scales = "free_y") +
#   theme_minimal() +
#   labs(x = "Station Code", y = "Bacteria Lysed", fill = "Analysis Method",
#        title = "Bacteria Lysed per Burst Size",
#        subtitle = "Faceted by VP Type and Burst Size") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# # Ensure the data has a 'bacteria_lysogeny' column calculated from lysogenic types
# vp_pe_df_combined <- vp_pe_df_combined %>%
#   mutate(bacteria_lysogeny = if_else(VP_type == "Lysogenic", percent_bacteria_lys, NA_real_))
# 
# # Prepare data for plotting
# plot_data <- vp_pe_df_combined %>%
#   filter(VP_type %in% c("Lytic", "Lysogenic")) %>%
#   mutate(station_code = factor(station_code, levels = unique(station_code)))
# 
# # Plotting
# ggplot(plot_data, aes(x = station_code, y = percent_bacteria_lys, fill = analysis_method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_grid(burst_size ~ VP_type, scales = "free_y") +
#   theme_minimal() +
#   labs(x = "Station Code", y = "Bacteria Lysed", fill = "Analysis Method",
#        title = "Bacteria Lysed per Burst Size",
#        subtitle = "Faceted by VP Type and Burst Size") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


