source("scripts/0_source.R")

filtered_counts <- read.csv( "results/viral_loss_corrected_counts/FILTERED_PE_Cruises_FCS_VIRAL_LOSS_CORRECTED_0_24_per_mL.csv")
bacterial_growth <- filtered_counts %>%
  group_by(Location_Station, Timepoint) %>%
  summarise(mean_c_Bacteria = mean(c_Bacteria, na.rm = TRUE)) %>%
  arrange(Location_Station, Timepoint)
bacterial_growth <- bacterial_growth %>%
  group_by(Location_Station)%>%
  mutate(
    growth_rate = (log(mean_c_Bacteria) - lag(log(mean_c_Bacteria))) / (Timepoint - lag(Timepoint)),
    generation_time = ifelse(!is.na(growth_rate) & growth_rate > 0, log(2) / growth_rate, NA)
  )
