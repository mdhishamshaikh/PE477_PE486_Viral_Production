# AIM: To create barplots for all variables ####

# 0.0 Setting up #####
source("scripts/0_source.R")

# 1.0 Importing combined PE data frame ####
pe_df <- read.csv( "./results/PE477_PE486_3depths_combined.csv") %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_")) %>%
  dplyr::filter(!Location_Station %in% c("PE477_7", "PE486_8")) # No VP assay was performed

station_variables <- c("Location", "Station_Number", "Depth", "sample_tag", "Location_Station")

# 2.0 Plotting variables that have information at 3 depths ####

variables_3_depths <- c("Salinity", "Temperature", "Density", "Conductivity", "Turbidity", "Oxygen", "Fluorescence",
                        "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
                        "Total_Viruses", "V1", "V2", "V3",
                        "VBR",
                        "Nitrate", "Nitrite", "Phosphate", "Silicate")


  
pe_df_long_3_depths <- pe_df %>%
  dplyr::select(all_of(c(station_variables, variables_3_depths))) %>%
  tidyr::pivot_longer(cols = -all_of(station_variables),
                      names_to = "Variable",
                      values_to = "Value") 




# Split the data by 'Variable'
plot_list <- split(pe_df_long_3_depths, pe_df_long_3_depths$Variable)

# Create one plot per variable, with Depth as a facet
plots <- lapply(names(plot_list), function(var) {
  ggplot(plot_list[[var]], aes(x = Location_Station, y = Value, fill = 'black')) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, fill = 'black') +
    labs(title = paste(var), x = "Station", y = "Value") +
    facet_grid( Depth ~ Location, scales = "fixed") +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    theme_test(base_size =12) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside")
})
plots[[5]]


output_folder <- "./figures/variables_3_depths"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Loop through the plots and save them
for (i in seq_along(plots)) {
  # Generate a filename for each plot
  plot_name <- paste0(output_folder, "/", names(plot_list)[i], ".png")
  
  # Save the plot
  ggsave(
    filename = plot_name,
    plot = plots[[i]],
    width = 4, height = 6, dpi = 300  # Adjust dimensions and resolution as needed
  )
}


# 3.0 Plots for 7 m depth ####


variables_7m <- c("Salinity", "Temperature", "Density", "Conductivity", "Turbidity", "Oxygen", "Fluorescence",
                        "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
                        "Total_Viruses", "V1", "V2", "V3",
                        "VBR",
                        "Nitrate", "Nitrite", "Phosphate", "Silicate",
                  "VP_Lytic", "VP_Lysogenic", "decay_rate_exponential")


pe_df_long_7m <- pe_df %>%
  dplyr::select(all_of(c(station_variables, variables_7m))) %>%
  tidyr::pivot_longer(cols = -all_of(station_variables),
                      names_to = "Variable",
                      values_to = "Value") 




# Split the data by 'Variable'
plot_list_7m <- split(pe_df_long_7m, pe_df_long_7m$Variable)

# Create one plot per variable, with Depth as a facet
plots_7m <- lapply(names(plot_list_7m), function(var) {
  ggplot(plot_list_7m[[var]], aes(x = as.factor(Station_Number), y = Value, fill = 'black')) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, fill = 'black') +
    labs(title = paste(var), x = "Station", y = "Value") +
    facet_grid( ~ Location, scales = "fixed") +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    theme_test(base_size =8) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside")
})
plots_7m[[5]]


output_folder <- "./figures/variables_7m"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Loop through the plots and save them
for (i in seq_along(plots_7m)) {
  # Generate a filename for each plot
  plot_name <- paste0(output_folder, "/", names(plot_list_7m)[i], ".png")
  
  # Save the plot
  ggsave(
    filename = plot_name,
    plot = plots_7m[[i]],
    width = 4, height = 3, dpi = 300  # Adjust dimensions and resolution as needed
  )
}

