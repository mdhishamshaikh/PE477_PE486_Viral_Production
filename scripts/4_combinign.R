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

# The sample_metadata, abundance, nutrients, salinity dataframe is ready 


#### 2.0 Net bacterial growth rate per station ####
net_growth <- read.csv("PE_Cruises/results/PE_Cruises_viral_production/vp_BP_calc.csv")

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


#### 3.0 Viral decay rates 

viral_decay <- read.csv("./results/viral_decay.csv")
