library(tidyverse)
library(devtools)
library(readxl)

#Install viral production calculator from Github
#install_github("mdhishamshaikh/ViralProduction_R")
library(viralprod)
#issues with lme4 package. It's important that we install it form source
#install.packages("lme4", type = "source") 
library(lme4)


setwd("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production")
#Import FCM count csv
counts_per_mL<- read.csv("./PE_Cruises/results/PE_Cruises_per_mL.csv")
str(counts_per_mL)
counts_per_mL<- counts_per_mL %>%
  dplyr::filter(counts_per_mL$Sample_Type != '0.22')

#There's some additional files that I would like to exclude from this. These are replicates with low events/sec
#To do so I will use the selected data from the metadata file. This should have already been done before processing. 
#Importing the .xlsx file

selected_files<- read_excel("Metadata/PE_Cruises_VP_Metadata.xlsx", sheet = 'Selected_Metadata') #contains selected filenames

counts_per_mL<- inner_join(counts_per_mL, selected_files, by = "Sample_Name")
#write.csv(counts_per_mL, file = 'Linear_Regression_Method_Supervised/PE477_PE486_filtered.csv', row.names = F)


 #Import abundance data
url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
str(abundance)

#Extracting first entries per Location/Station_Number combinations, as some stations had multiple depths and this information is missing from the fiel on GitHub.
#The first one is the depth for VP assays.
abundance<- abundance %>%
  filter(Depth %in% c(1, 7) )

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
              methods = c(2,10),
              write_output = T,
              output_dir = 'PE_Cruises2/results/PE_Cruises_viral_production')

# 
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


