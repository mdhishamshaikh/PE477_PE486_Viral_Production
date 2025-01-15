# AIM: Importing microbial (bacterial, viral and phytoplankton) abundance and nutrients data from `github.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/`

# 0.0 Setting up ####
source("./scripts/0_source.R")

# 1.0 Importing bacterial and viral abundances along with nutrients
bv_nuts <- readr::read_csv("https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv") %>%
  dplyr::filter(Location != "NJ2020")

# 2.0 Importing phytoplankton abundance
phyto <- readr::read_csv("https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/phytoplankton_counts_pe477_pe486.csv")
# phyto counts are not in per mL right now

# extracting cyanobacterial gate 14

cyano <- phyto %>%
  dplyr::filter(Gate == 14) %>%
  dplyr::select(c("Total_Events", "Location", "Station_Number", "Depth")) %>%
  dplyr::rename(Cyanobacteria = Total_Events)

# 3.0 Combining the two data frames and saving it ####

microbial_abundance_nutrients <- bv_nuts %>%
  dplyr::left_join(cyano, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::select(-c("Temperature", "Salinity", "TON", ends_with("Sample_Name")))

# Saving
dir.create("./results/microbial_abundance_nutrients/", recursive = T)
write.csv(microbial_abundance_nutrients, "./results/microbial_abundance_nutrients/microbial_abundance_nutrients.csv", row.names = F)
