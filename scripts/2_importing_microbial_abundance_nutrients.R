# AIM: Importing microbial (bacterial, viral and phytoplankton) abundance and nutrients data from `github.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/`

# 0.0 Setting up ####
source("./scripts/0_source.R")

# 1.0 Importing bacterial and viral abundances along with nutrients
bv_nuts <- readr::read_csv("https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv") %>%
  dplyr::filter(Location != "NJ2020")

# 2.0 Importing phytoplankton abundance
phyto <- readr::read_csv("https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/phytoplankton_counts_pe477_pe486.csv")



# extracting total phytoplankton gate 20

total_phyto <- phyto %>%
  dplyr::filter(Gate == 20) %>%
  dplyr::select(c("cells_per_mL", "Location", "Station_Number", "Depth")) %>%
  dplyr::rename(Total_phyto = cells_per_mL) %>%
  distinct()


# extracting cyanobacterial gate 14

cyano <- phyto %>%
  dplyr::filter(Gate == 14) %>%
  dplyr::select(c("cells_per_mL", "Location", "Station_Number", "Depth")) %>%
  dplyr::rename(Cyanobacteria = cells_per_mL) %>%
  distinct()


# extracting doggers phyto 1 gate 18

doggers_phyto_1 <- phyto %>%
  dplyr::filter(Gate == 18) %>%
  dplyr::select(c("cells_per_mL", "Location", "Station_Number", "Depth")) %>%
  dplyr::rename(Doggers_phyto_1 = cells_per_mL) %>%
  distinct()

# extracting doggers phyto 2 gate 19

doggers_phyto_2 <- phyto %>%
  dplyr::filter(Gate == 19) %>%
  dplyr::select(c("cells_per_mL", "Location", "Station_Number", "Depth")) %>%
  dplyr::rename(Doggers_phyto_2 = cells_per_mL) %>%
  distinct()


# 3.0 Combining the two data frames and saving it ####

# replacing zeros with half of detecting limit
dl_nitrate <- 0.1529
dl_phosphate <- 0.0152
dl_silicate <- 0.0711

microbial_abundance_nutrients <- bv_nuts %>%
  dplyr::left_join(total_phyto, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::left_join(cyano, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::left_join(doggers_phyto_1, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::left_join(doggers_phyto_2, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::select(-c("Temperature", "Salinity", "TON", ends_with("Sample_Name"))) %>%
  mutate(
    Nitrate = ifelse(Nitrate == 0, dl_nitrate/2, Nitrate),
    Silicate = ifelse(Silicate == 0, dl_phosphate/2, Silicate),
    Phosphate = ifelse(Phosphate == 0, dl_silicate/2, Phosphate),
    Nitrite = ifelse(Nitrite < 0, 0, Nitrite)
  )



# Saving
dir.create("./results/microbial_abundance_nutrients/", recursive = T)
write.csv(microbial_abundance_nutrients, "./results/microbial_abundance_nutrients/microbial_abundance_nutrients.csv", row.names = F)
