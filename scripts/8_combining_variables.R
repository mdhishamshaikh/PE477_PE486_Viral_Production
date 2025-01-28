# AIM: To combine all the variables in one dataframe ##

# 0.0 Setting up ####
source("scripts/0_source.R")


# 1.0 CTD profiles #
# Here I will only work with CTD variables from 3 depths (7 m, 15 m and 30 m)
# Variables: Salinity, temperature, pressure, density, conductivity, oxyge, turbidity, fluorescence
ctd <- read.csv("./results/ctd_profiles/ctd_variables_sampling_depths.csv") %>%
  dplyr::select(-c(Measured_depth, Distance, Location_Station_Number)) %>%
  dplyr::rename(Depth = Target_Depth) 


# 2.0 Microbial abundances and nutrients ####

abundance_nuts <- read.csv("./results/microbial_abundance_nutrients/microbial_abundance_nutrients.csv")
# Missing 15 m and 30 m from PE477_1 and PE477_2 as samples were not taken for them.


# 3.0 Viral production ####

vp <- read.csv("./results/viral_production_analyses/viralproduction_viralprod_lr.csv") %>%
  dplyr::filter(VP_Method == "VIPCAL-SE") %>% # only retaining VIPCAL-SE
  tidyr::pivot_wider(names_from = Treatment,
                     values_from = c(VP, abs_VP, VP_SE),
                     names_prefix = "") %>%
  dplyr::select(-c(Time_Range, Population, Location_Station, VP_Method)) 
#As Time_Range is T0_T24 for all, population we focus on is c_Viruses, and VIPCAL-SE is the only method
# VP_Lytic and VP_Lysogenic still needs to be corrected for bacterial loss when processing
#VP_Lytic also needs to be corrected for viral loss due t sticking on the wall

# 4.0 Viral decay ####

decay <- read.csv("./results/viral_decay.csv") %>%
  dplyr::filter(Population == "c_Viruses") %>%
  dplyr::select(-c(ends_with("exponential"), Sample_Type, Population))

# 5.0 Correcting lytic and lysogenic viral production for bacterial loss & loss from viral decay ####

vp_corrected_for_decay <- vp %>%
  left_join(decay, by = c("Location", "Station_Number", "Depth")) %>%
  mutate(Corrected_VP_Lytic = (VP_Lytic + decay_rate_linear)*100/bac_efficiency, # Also correcting viral loss for bacterial loss during processing
         Corrected_VP_Lysogenic = (VP_Lysogenic)*100/bac_efficiency) # Correctin gonly for bacterial loss


# 5.0 Combining all variables ####

pe_df <-  ctd %>%
  left_join(abundance_nuts, by = c("Location", "Station_Number", "Depth")) %>%
  left_join(vp_corrected_for_decay, by = c("Location", "Station_Number", "Depth")) %>%
  mutate(sample_tag = paste(Location, Station_Number, Depth, sep = "_"))

# 6.0 Calculating bacterial loss and percent lysogeny ####

burst_size <- 17, 40, 50
pe_df <- pe_df %>%
  mutate(bacterial_loss_rate = Corrected_VP_Lytic/burst_size, # assuming burst size of 25
         percent_bacterial_loss_day = (bacterial_loss_rate*24*100)/Total_Bacteria, # The percent of standing stock lysed over a day
         percent_lysogeny = (Corrected_VP_Lysogenic*100)/(burst_size*Total_Bacteria), 
         viral_turnover_time = Total_Viruses/Corrected_VP_Lytic # Time taken to replace the current standing stock of viruses
  )
plot(pe_df$percent_bacterial_loss_day)
plot(pe_df$percent_lysogeny)

burst_sizes <- c(17, 40, 50)

# Loop through burst sizes and dynamically add columns
for (b in burst_sizes) {
  pe_df <- pe_df %>%
    mutate(
      !!paste0("bacterial_loss_rate_burst_", b) := Corrected_VP_Lytic / b,
      !!paste0("percent_bacterial_loss_day_burst_", b) := (Corrected_VP_Lytic / b * 24 * 100) / Total_Bacteria,
      !!paste0("percent_lysogeny_burst_", b) := (Corrected_VP_Lysogenic * 100) / (b * Total_Bacteria)
    )
}

# Plot one of the new dynamic columns as an example
plot(pe_df$bacterial_loss_rate_burst_17)
plot(pe_df$percent_bacterial_loss_day_burst_40)


write.csv(pe_df, "./results/PE477_PE486_3depths_combined.csv", row.names = F)

pe_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7)

write.csv(pe_df_7m, "./results/PE477_PE486_3depths_combined_7m.csv", row.names = F)



ggplot(pe_df %>% dplyr::filter(Depth == 7), aes(x = as.factor(sample_tag), y = Corrected_VP_Lytic , fill = 'black')) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, fill = 'black') +
  #labs(title = paste(var), x = "Station", y = "Value") +
  #facet_grid( ~ Location, scales = "fixed") +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  theme_test(base_size =14) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside")

