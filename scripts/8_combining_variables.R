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


# 4.0 Viral decay ####

decay <- read.csv("./results/viral_decay.csv") %>%
  dplyr::filter(Population == "c_Viruses") %>%
  dplyr::select(-c(ends_with("linear"), Sample_Type, Population))


# 5.0 Combining all variables ####

pe_df <-  ctd %>%
  left_join(abundance_nuts, by = c("Location", "Station_Number", "Depth")) %>%
  left_join(vp, by = c("Location", "Station_Number", "Depth")) %>%
  left_join(decay, by = c("Location", "Station_Number", "Depth")) %>%
  mutate(sample_tag = paste(Location, Station_Number, Depth, sep = "_"))

write.csv(pe_df, "./results/PE477_PE486_3depths_combined.csv", row.names = F)
