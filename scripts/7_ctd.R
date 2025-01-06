# The aim of the script is to extract clusters from the CTD profiles. 

# 0.0 Setting up ####
library(tidyverse)
library(ggsci)
library(vegan)
library(ggfortify)
library(psych)

# 1.0 Importing CTD profiles.
ctd_profiles <- utils::read.csv("./results/ctd_profiles/ctd_profiles.csv")

# 2.0 Factor analysis ####

# KMO and Barlett's test to check if one can perform fcator analysis
kmo_result <- KMO(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, -Pressure))
bartlett_result <- cortest.bartlett(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, - Pressure))

print(kmo_result)
print(bartlett_result)
# The data is suitable. But as Salinity has low MSA, i will exclude it. 

kmo_result_no_salinity <- KMO(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, -Pressure, -Salinity))
bartlett_result_no_salinity <- cortest.bartlett(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, - Pressure, -Salinity))

print(kmo_result_no_salinity)
print(bartlett_result_no_salinity)
# Improves the values a little bit.

# Factor analysis
fa_result <- fa(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth, -Pressure, -Salinity), 
                nfactors = 2, rotate = "varimax", fm = "ml")
print(fa_result)


