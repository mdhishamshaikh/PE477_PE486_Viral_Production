library(tidyverse)
library(ggsci)
library(readxl)
library(viralprod)

cruise_vp<- read.csv("PE_Cruises/results/PE_Cruises_viral_production/vp_results_BP.csv")


# I will work on total viruses only for now. I will combine this with the viral production rates calculated using supervised linear regression model.

vp_df <- cruise_vp %>%
  dplyr::filter(!Sample_Type %in% c("0.22", "VPC"),
                Population == 'c_Viruses',
                VP_Method == 'VPCL_AR_DIFF_SE')  %>%
  mutate(Depth = 7)


# importing LM supervised data

lm_supervised <- read_excel("Linear_Regression_Method_Supervised/PE477_PE486_VP_LM.xlsx", 
                            sheet = "final_output") %>%
  select(!c(station_code, corrected_Lytic, corrected_Lysogenic, bac_efficiency)) %>%
  rename(VP = Lytic, Diff = Lysogenic) %>%
  pivot_longer(cols = c(VP, Diff), names_to = "Sample_Type", values_to = "VP") %>%
  mutate(VP_Method = "LM_Supervised",
         Time_Range = NA,
         Population = "c_Viruses",
         abs_VP = NA,
         VP_SE = NA,
         VP_R_Squared = NA)
 
  

vp_lm_df <- rbind(vp_df, lm_supervised)  %>%
  mutate(Depth = 1)


# Importing abudnance data

url<- "https://raw.githubusercontent.com/mdhishamshaikh/NorthSea_Microbial_Abundance_Nutrients/main/nj2020_pe477_pe486_bv_abundance_abiotic.csv"
abundance<- readr::read_csv(url)
abundance<- abundance %>%
  filter(Depth %in% c(1, 7),
         Location != "NJ2020") %>%
  mutate(Depth = 1)  

#Some checks
vp_check_populations(counts_per_mL)

class(abundance)
abundance <- vp_class_ori_abu(abundance)
class(abundance)

class(counts_per_mL) #failed
counts_per_mL <- vp_class_count_data(counts_per_mL) #passed
class(counts_per_mL)
# Running viral prod vp_analyze

viralprod::vp_analyze(x = counts_per_mL,
                      vp_results = vp_lm_df,
                      output_dir = "./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM",
                      original_abundances = abundance)



# Importing vp_analyzed file

vp_analyzed <- read.csv("./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_results_ANALYZED.csv")

vp_analyzed <- vp_analyzed %>% 
  mutate(
    Sample_Type = recode(Sample_Type, "VP" = "Lytic", "Diff" = "Lysogenic"),
    VP_Method = recode(VP_Method, "VPCL_AR_DIFF_SE" = "VIPCAL-SE", "LM_Supervised" = "LM-S")
  )
  
# Recode VP_Method
vp_analyzed <- vp_analyzed %>%
  mutate(VP_Method = recode(VP_Method, "VIPCAL-SE" = "VPCL", "LM-S" = "LM")) %>%
  select(-c("Time_Range", "VP", "abs_VP", "VP_SE", "VP_R_Squared")) 


# Columns to keep without prefix
keep_cols <- c("Location", "Station_Number", "Depth", "Population", "Sample_Type", "VP_Method", "B_0", "B_OS", "V_OS" )

# Columns to pivot
pivot_cols <- setdiff(colnames(vp_analyzed), keep_cols)

# Pivot wider with prefixes from Sample_Type and VP_Method
vp_analyzed_wide <- vp_analyzed %>%
  pivot_longer(cols = all_of(pivot_cols), names_to = "Variable", values_to = "Value") %>%
  unite("Prefix", Sample_Type, VP_Method, Variable, sep = "_") %>%
  pivot_wider(names_from = Prefix, values_from = Value)
vp_analyzed_wide_ts <- vp_analyzed_wide %>%
  left_join(abundance, by = c("Location", "Station_Number", "Depth")) %>%
  mutate(Location2 = Location,
         Station_Number2 = Station_Number)


write.table(vp_analyzed_wide_ts, "./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_analyzed_pe_cruises_odv_input.txt", sep = "\t",  row.names = F, quote = F)
