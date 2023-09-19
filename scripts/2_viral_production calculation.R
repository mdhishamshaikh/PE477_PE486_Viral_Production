library(tidyverse)
library(devtools)

#Install viral production calculator from Github
install_github("mdhishamshaikh/ViralProduction_R")
library(viralprod)

#Import FCM count csv

counts_per_mL<- read.csv("./PE_Cruises/results/PE_Cruises_per_mL.csv")
str(counts_per_mL)
counts_per_mL<- counts_per_mL %>%
  dplyr::filter(counts_per_mL$Sample_Type != '0.22')

#Import abundance data
url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
str(abundance)  



#Some checks
vp_check_populations(counts_per_mL)
vp_class_ori_abu(abundance) #failed
#making necessary changes
colnames(abundance)[2]<- "Station_Number"
abundance$Station_Number <- as.integer(abundance$Station_Number)
vp_class_ori_abu(abundance) #passed


class(abundance)
class(vp_class_ori_abu(abundance))

vp_class_count_data(counts_per_mL) #passed


#Running viralprod

vp_end_to_end(data = counts_per_mL,
              original_abundances = abundance,
              write_output = T,
              output_dir = 'PE_Cruises/results/PE_Cruises_viral_production')

#vp_visualize failed

#### Visualize

cruise_vp<- read.csv("PE_Cruises/results/PE_Cruises_viral_production/vp_results_BP.csv")
cruise_vp<- cruise_vp %>%
  dplyr::filter(VP_Method == 'VPCL_AR_DIFF_LMER_SE',
                Sample_Type != 'VPC',
                Population == 'c_Viruses') %>%
  mutate(Stations = paste0(Location, "_",as.character(Station_Number)))
str(cruise_vp)

cruise_vp

ggplot(data = cruise_vp,
       aes(x = Stations,
           y = abs_VP,
           color = Sample_Type,
           fill = Time_Range))+
  geom_point()


abundance7<- abundance %>%
  dplyr::filter(Depth %in% c(1,7),
                Location != "NJ2020") %>%
  mutate(Depth = as.integer(Depth))
str(abundance7)

vp<- vp_results_output_BP_df %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_LMER_SE",
                Sample_Type != 'VPC',
                Population == 'c_Viruses') %>%
  mutate(Sample_Type = base::ifelse(Sample_Type == 'VP', 'Lytic', 'Lysogenic')) %>%
  select(-c(VP_Method, VP_R_Squared, VP_SE)) %>%
  pivot_wider(names_from = Sample_Type,
              values_from = c(VP, abs_VP))
vp$Depth<- 7
str(vp)
#need to wrangle it tolytic andlysogenic prod columns

vp_abundance <- full_join(vp, abundance7, by = c("Location", "Station_Number", "Depth"))
head(vp_abundance)

#REMOVE ROWS THAT DON'T HAVE VP assay perfomred.

vp_abundance<- vp_abundance %>% 
  dplyr::filter(!(is.na(VP_Lytic)))

#This is the combined datframe. Save it as csv and we'll use thisf or visualization.

write.csv(vp_abundance, paste0(getwd(), 
                               "/PE_Cruises/results/PE_Cruises_viral_production/vp_abundance_nutrients_ts.csv"),
          row.names = F)
