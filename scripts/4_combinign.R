# The aim fo this script is to combine all the dataI have

# 0.0 Setting up ##
library(tidyverse)
library(ggsci)
library(lubridate)

#### 1.0 Sample metadata, Abundance, Nutrients, Temperature, and Salinity ####

# Importing metadata file
abundance_ts <- read.csv("./Data/nutrients_ts_abundance_coordinates/nutrients_ts_abundance_coordinates.csv")

# Removing bacterial and viral sample names and TON
abundance_ts <- abundance_ts %>%
  select(!ends_with("Sample_Name"), -TON)

# Formatting the date column
abundance_ts$Expt_Date <-lubridate::dmy(abundance_ts$Expt_Date)



# Replacing negative nutrients values with zero

columns_to_modify <- c("Nitrate", "Nitrite", "Silicate", "Phosphate")
abundance_ts[columns_to_modify] <- lapply(abundance_ts[columns_to_modify], function(x) ifelse(x < 0, 0, x))



# The sample_metadata, abundance, nutrients, salinity dataframe is ready 


#### 2.0 Net bacterial growth rate per station ####
net_growth <- read.csv("PE_Cruises_viral_production_OUTLIERS/vp_BP_calc.csv")

net_growth<- net_growth %>%
  group_by(combi_tag) %>%
  filter(Timepoint == max(Timepoint)) %>% # max time point identified for that station
  ungroup() 
{ #quick plot
  
  net_growth_long <- net_growth %>%
    pivot_longer(cols = c(Generation_Time, Net_Growth_Rate), names_to = "Metric", values_to = "Value")
  
  ggplot(net_growth_long, aes(x = combi_tag, y = Value, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    labs(x = "Combi Tag", y = "Value") +
    theme_minimal() +
    scale_fill_npg() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Generation Time and Net Growth Rate by Combi Tag and Population")  
  }

net_growth_wide <- net_growth %>%
  pivot_wider(
    names_from = Population,
    values_from = c(Generation_Time, Net_Growth_Rate)
  ) %>%
  select(-combi_tag) %>%
  select(-Timepoint) # also getting rid of timepoint as it is the same in all, and should be the same for comparing net bacterial growth rates

# Higher growth was observed in HNA as compared to LNA in the assay. 
# Bacterial growth dataframe is ready 

vp_pe_df <- abundance_ts %>%
  left_join(net_growth_wide, by = c("Location", "Station_Number", "Depth"))

#### 3.0 Viral decay rates 

viral_decay <- read.csv("./results/viral_decay.csv")

viral_decay_wide <- viral_decay %>%
  pivot_wider(
    names_from = Population,
    values_from = c(decay_rate, percent_decay_day)
  ) %>%
  select(-Sample_Type)

vp_pe_df <- vp_pe_df %>%
  left_join(viral_decay_wide, by = c("Location", "Station_Number", "Depth"))

#### 4.0 Viral Production data ####

# vp <- read.csv("./PE_Cruises_viral_production_OUTLIERS/vp_results_BP.csv", sep = ",")
# 
# vp <- vp %>%
#   select(-contains("BP")) %>% # Removing bacterial production data, as we did not provide any. I could also provide b
# select(-contains("abs")) # removing absolute viral production
# 
# 
# # Removing all the nutrient release from lysogenic induction calculations
# colnames_vp <- colnames(vp)
# 
# vp <- vp %>%
# select(-all_of(colnames_vp[grep("Lysogenic.*(DOC|DON|DOP|V_TT)", colnames_vp)]))
# 
# colnames_vp <- colnames(vp)
# colnames_vp <- gsub("P_B_Loss", "Percent_Bacteria_Loss", colnames_vp)
# colnames_vp <- gsub("P_Cells", "Percent_Cells", colnames_vp)
# colnames_vp <- gsub("Rate", "Bacteria_Loss_Rate", colnames_vp)
# colnames_vp <- gsub("TT", "Turnover_Time", colnames_vp)
# colnames(vp) <- colnames_vp
vp <- read.csv("./results/viral_production_vipcal_se.csv", sep = ",")
vp <- vp %>%
  mutate(Depth = 1)

vp_pe_df <- vp_pe_df %>%
  left_join(vp, by = c("Location", "Station_Number", "Depth"))


#### 5.0 Adding grazing due to heteronanoflagellate data ####

flb <- read.csv("./results/flb_grazing.csv")
# Removing vcontrols for now

flb <- flb %>%
  filter(Sample_Type != "Control") %>%
  select(-Sample_Type)

vp_pe_df <- vp_pe_df %>%
  left_join(flb, by = c("Location", "Station_Number", "Depth"))



#### 6.0 Saving it for ODV ####

vp_pe_df_odv <- vp_pe_df %>%
  mutate(Location2 = Location,
       Station_Number2 = Station_Number,
       HNALNA = HNA/LNA,
       V1V2 = V1/V2,
       V1V3 = V1/V3,
       V2V3 = V2/V3) %>%
  filter(!is.na(Generation_Time_c_Bacteria)) %>%# remove abuudance and TS data for the extra stations where VP assay wasn't performed.
dplyr::rename(Cruise = Location,
       Station = Station_Number
       ) %>%
  select(Cruise, Station, everything())

write.table(vp_pe_df_odv, "./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_analyzed_pe_cruises_combined.txt", sep = "\t",  row.names = F, quote = F)



