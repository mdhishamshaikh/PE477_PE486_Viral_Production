
# Clcualting burs sizes

counts_per_mL<- read.csv("./PE_Cruises_FCS/results/PE_Cruises_FCS_per_mL.csv")
str(counts_per_mL)

vpc <- counts_per_mL %>%
  dplyr::filter(Sample_Type == 'VPC',
                Staining_Protocol == "Viruses") %>%
  dplyr::select(c(Location, Station_Number, Timepoint, Replicate, c_Bacteria, c_Viruses)) %>%
  group_by(Location, Station_Number, Timepoint, Replicate) %>%
  slice_head(n = 1) %>%  # Retain the first row of each group
  ungroup()



burst_size <- vpc %>%
  group_by(Location, Station_Number, Replicate) %>%
  arrange(Timepoint) %>%  # Ensure timepoints are in order
  mutate(
    ΔViruses = c_Viruses - lag(c_Viruses),       # Change in viruses
    ΔBacteria = lag(c_Bacteria) - c_Bacteria,   # Change in bacteria
    Burst_Size = ΔViruses / ΔBacteria           # Burst size
  ) %>%
  ungroup()



library(dplyr)

# Calculate burst size between 0-hour and 24-hour values
burst_size_24hr <- vpc %>%
  dplyr::filter(Timepoint %in% c(0, 24)) %>%  # Retain only 0-hour and 24-hour timepoints
  arrange(Location, Station_Number, Replicate, Timepoint) %>%  # Ensure correct order
  group_by(Location, Station_Number, Replicate) %>%
  summarise(
    ΔViruses = c_Viruses[Timepoint == 24] - c_Viruses[Timepoint == 0],  # Change in viruses
    ΔBacteria = c_Bacteria[Timepoint == 0] - c_Bacteria[Timepoint == 24],  # Change in bacteria
    Burst_Size = ΔViruses / ΔBacteria,  # Burst size calculation
    .groups = "drop"
  )

# Check the result
head(burst_size_24hr)


#FLUCTUATIONS AND LACK OF AGREEMENT BETWEEN REPLCIATES SUGGEST THIS IS NOT THE BEST WAY TO CALCULATE BURST SIZE