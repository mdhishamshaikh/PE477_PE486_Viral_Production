# Setting up ####

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv")

# 2.0 PCA - physicochemical parameters ####

physicochemical_params <- c(Temperature, Salinity, Density, Conductivity, Turbidity,
                            Nitrate, Phosphate, Silicate)
bio_params <- c(Oxygen, Fluorescence,
                Total_Bacteria, HNA, LNA, Cyanobacteria,
                Total_Viruses, V1, V2, V3,
                VBR,
                )
sample_params <- c(Location, Station_Number, Depth, 
                   sample_tag, Latitude, Longiude)