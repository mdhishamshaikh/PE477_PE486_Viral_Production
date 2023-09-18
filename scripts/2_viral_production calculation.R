library(tidyverse)
library(devtools)

#Install viral production calculator from Github
install_github("mdhishamshaikh/ViralProduction_R")
library(viralprod)

#Import FCM count csv

counts_per_mL<- read.csv("./PE_Cruises/results/PE_Cruises_per_mL.csv")
str(counts_per_mL)

#Import abundance data
