{
library(tidyverse)
library(ggsci)
library(vegan)
library(ggfortify)
  library(psych)
}


ctd_profiles <- utils::read.csv("./results/ctd_profiles/ctd_profiles.csv")

# 1.0 PE477 & PE486 combined ####

# 1.1 Factor analysis on the entire dataset. ####

# Calculating KMO and Bartlett's Test to check suitability of data
kmo_result <- KMO(log(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth)))
bartlett_result <- cortest.bartlett(ctd_profiles %>% select(-Location_Station_Number, -Location, -Station_Number, -Depth))

print(kmo_result)
The Kaiser

print(bartlett_result)

