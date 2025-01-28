
#The aim of this script is to process `fcs` files and metadata to extract per mL values for PE477 and PE486 cruises

# Steps include:\
# 1. Defining general variables.\
# 2. Setting up the environment.\
# 3. Importing and processing metadata file.\
# 4. Importing FCS files.\
# 4. Creating a reference FCS file.\
# 5. Perform a small test.\
# 6. Extract cytoplots.\
# 7. Extract counts.\
# 8. Ajust TE and tidy data.\
# 
# 

#Step 1: Setting up ####
source("scripts/0_source.R")

#Define the project title
project_title<- "PE_Cruises_FCS_0.22"


# Define the working directory `work_dir`.\
# This could be done by using `getwd()` function.\

proj_dir <- getwd() # Making sure it is correct


# **2. Setting up environment**\ ####
# This is to create a space for your data, metadata, raw data
# along with a directory for output.\



#Setting up the directories using a project name
set_up_vp_count(project_title)
# Also sets up the `work_dir`. This is the directory where everything related to fcs counting will be stored.

# **3.  Import and process metadata file**\ ####
# The function `metadata_processing` accepts file path to the metadata file. It 
# also takes care of the dates.
# Please ensure that you have the following parameters mentioned in your metadata file.\
# ADD THE PARAMETERS HERE\

metadata_processing(paste0(proj_dir,"/metadata/PE_Cruises_VP_Metadata.xlsx"), sheet = "Selected_Metadata", extension = ".xlsx", project_title = "PE_Cruises")

# #I will susbet metadata for 0.22 and TE
metadata<- metadata %>%
  dplyr::filter(Sample_Type %in% c("0.22", "TE"),
                Staining_Protocol == "Viruses")
#metadata <- metadata[1:10,]

# **4. Importing FCS files**\ ####
# We can utilize the `Sample_Name` variable in the `metadata` file to copy files for analyses.\
# We need to provide the origin of these FCS files. Files will be copied to `./data/raw_data` under proj_dir\
# 

# Data is stored under `data` folder for now. This can also have file swe do not need.
fcs_origin<- paste0(proj_dir,"/data/") # Use absolute path

# The `import_fcs` functions accepts fcs file origin directory project title (already defined),  and moves the files to `proj_dir`/data/raw_data
import_fcs(fcs_origin)



#### 5.0 Gating ####

#Let's populate the gating dataframe with default values.
#The default gating values

{
  metadata = metadata
  b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6)
  b_fl1 = c(1.9, 1.4, 1.4,1.8, 2.8, 3.7, 3.2, 2.75)
  v_ssc = c(-0.25, 1.2)
  v_fl1 = c(-0.1, 1.7)
  hna_ssc = c(0.6, 3.5)
  hna_fl1 = c(2.15, 3.5)
  lna_ssc = c(1.0, 3.5)
  lna_fl1 = c(1.4,2.15)  
  v1_ssc = c(-0.1, 0.90)
  v1_fl1 = c(-0.1, 0.8)
  v2_ssc = c(-0.1, 0.90)
  v2_fl1 = c(0.8, 1.25)
  v3_ssc = c(-0.1, 1.3)
  v3_fl1 = c(1.25, 1.7)
}


populate_gate_df() #There are default gates for Bacteria (total, HNA, LNA) and viruses (total, V1, V2, V3)
#For now we use te values defined above
#Gates can be adjusted later if needed using the same function
#Now let's evaluate how it performs on our samples

#We always only start with plots

get_bv_plots()
# saved under results
#Note changes that need to be made. You can visuaize a cytoplot using sample index
#and change gate_df using populate_df to make changes. 



##### i will make these changes between V1 and V2 later. Fo rnow, total viruses are alright. 
#Let's get the stats out

get_bv_stats()


#Import csv
counts<- read.csv(paste0(project_title, "/results/", project_title, "_counts.csv"))
str(counts)


combine_metadata_counts(metadata_df = metadata, counts_df = counts) # Adds metadata
calc_TE(df = counts_metadata) # calculates TE for sample
adj_counts <- adjust_TE(counts_metadata_df = counts_metadata_TE) # saves per mL counts

# Importing output per mL file
counts_per_mL<- read.csv(paste0(project_title, "/results/", project_title, "_per_mL.csv"))

counts_per_mL_0.22 <- counts_per_mL %>%
  dplyr::filter(Sample_Type == "0.22")


# Create a combined variable for Location and Station_Number
counts_per_mL$Location_Station <- paste(counts_per_mL$Location, counts_per_mL$Station_Number, sep = "_")

counts_per_mL_outliers_removed <- counts_per_mL %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE477" & Station_Number == 2 & Timepoint == 6 & Replicate == 1,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  ))

counts_per_mL_outliers_removed$Location_Station <- paste(counts_per_mL_outliers_removed$Location, counts_per_mL_outliers_removed$Station_Number, sep = "_")


# Create the plot
ggplot(counts_per_mL_outliers_removed, aes(x = Timepoint, y = c_Viruses, color = as.factor(Replicate))) +
 # geom_line(aes(group = interaction(Location_Station, Replicate)), size = 1) +
  geom_smooth(method = 'lm', se= F) +
  geom_point(size = 3) +
  facet_wrap(~ Location_Station, scales = "free") +
  labs(
    title = "0.22 µm Viral Counts Over Time",
    x = "Timepoint",
    y = "Viral Counts (c_Viruses)",
    color = "Replicate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )

# Re running TE calcualtions with only samples processed by Kirsten
counts_metadata<- counts_metadata %>%
  dplyr::filter(Date_Measurement %in% c("2022-03-31", "2022-04-06", "2022-04-08", "2022-04-13"))


calc_TE(df = counts_metadata) # calculates TE for sample
adj_counts <- adjust_TE(counts_metadata_df = counts_metadata_TE) # saves per mL counts
# Importing output per mL file
counts_per_mL_KK<- read.csv(paste0(project_title, "/results/", project_title, "_per_mL.csv"))


dim(counts_per_mL_all)
dim(counts_per_mL_KK)

#Now comparin counts per Ml between all timepoints and K timepoints

counts_all_3_12 <- counts_per_mL_all %>%
  dplyr::filter(Timepoint %in% c(3, 6, 9, 12))
dim(counts_all_3_12)
dim(counts_per_mL_KK)

plot(counts_all_3_12$c_Viruses, counts_per_mL_KK$c_Viruses)
lm(counts_all_3_12$c_Viruses ~ counts_per_mL_KK$c_Viruses)

# No differenc ebetween KK only and combiign everything when adjusting for TE

# Now checking if slop eo fT0 & T24 is ame as T3-T12


# Step 1: Create a Location_Station_Number variable
counts_per_mL_all <- counts_per_mL_all %>%
  mutate(Location_Station_Number = paste(Location, Station_Number, sep = "_"))

# Step 2: Average `c_Viruses` by Location_Station_Number and Timepoint
averaged_data <- counts_per_mL_all %>%
  group_by(Location_Station_Number, Timepoint) %>%
  summarize(mean_c_Viruses = mean(c_Viruses, na.rm = TRUE), .groups = "drop")

# Step 3: Calculate slopes for each scenario
calculate_slope <- function(data, timepoints) {
  filtered_data <- data %>% dplyr::filter(Timepoint %in% timepoints)
  
  slopes <- filtered_data %>%
    group_by(Location_Station_Number) %>%
    summarize(
      slope = if(n() > 1) coef(lm(mean_c_Viruses ~ Timepoint))[2] else NA,
      .groups = "drop"
    )
  return(slopes)
}

# 3.1: Slope using all timepoints
slope_all <- calculate_slope(averaged_data, unique(averaged_data$Timepoint))

# 3.2: Slope using only t0 and t24
slope_t0_t24 <- calculate_slope(averaged_data, c(0, 24))

# 3.3: Slope using t3, t6, t9, and t12
slope_t3_t12 <- calculate_slope(averaged_data, c(3, 6, 9, 12))

# View results
print("Slopes using all timepoints")
print(slope_all)

print("Slopes using only t0 and t24")
print(slope_t0_t24)

print("Slopes using t3, t6, t9, and t12")
print(slope_t3_t12)



plot((slope_t0_t24 %>% dplyr::filter(!Location_Station_Number %in% c("PE477_3", "PE477_5", "PE477_6")))$slope, slope_t3_t12$slope)
lm((slope_t0_t24 %>% dplyr::filter(!Location_Station_Number %in% c("PE477_3", "PE477_5", "PE477_6")))$slope ~ slope_t3_t12$slope)



##Bacteria #####

ggplot(counts_per_mL_outliers_removed %>% dplyr::filter(Staining_Protocol == "Viruses"), aes(x = Timepoint, y = c_Bacteria, color = as.factor(Replicate))) +
  # geom_line(aes(group = interaction(Location_Station, Replicate)), size = 1) +
  geom_smooth(method = 'lm', se= F) +
  geom_point(size = 3) +
  facet_wrap(~ Location_Station, scales = "free") +
  labs(
    title = "0.22 µm Bacterial Counts Over Time",
    x = "Timepoint",
    y = "Bacterial Counts ",
    color = "Replicate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )



counts_per_mL_outliers_removed <- counts_per_mL_outliers_removed %>%
  left_join(
    abundance_nuts %>%
      select(Location, Station_Number, Total_Bacteria, Total_Viruses, HNA, LNA, V1, V2, V3),
    by = c("Location", "Station_Number")
  ) 



counts_per_mL_outliers_removed<- counts_per_mL_outliers_removed %>%
  dplyr::mutate(vir_eff_0.22 = c_Viruses*100/Total_Viruses,
                bac_eff_0.22 = c_Bacteria*100/Total_Bacteria,
                hna_eff_0.22 = c_HNA*100/HNA,
                lna_eff_0.22 = c_LNA*100/LNA,
                v1_eff_0.22 = c_V1*100/V1,
                v2_eff_0.22 = c_V2*100/V2,
                v3_eff_0.22 = c_V3*100/V3)


# Avergaing replciates
average_efficiency <- counts_per_mL_outliers_removed %>%
  dplyr::filter(Staining_Protocol == "Viruses") %>%
  group_by(Location, Station_Number, Timepoint) %>%
  summarise(
    avg_bac_eff = mean(bac_eff_0.22, na.rm = TRUE),
    avg_vir_eff = mean(vir_eff_0.22, na.rm = TRUE),
    avg_hna_eff = mean(hna_eff_0.22, na.rm = TRUE),
    avg_lna_eff = mean(lna_eff_0.22, na.rm = TRUE),
    avg_v1_eff = mean(v1_eff_0.22, na.rm = TRUE),
    avg_v2_eff = mean(v2_eff_0.22, na.rm = TRUE),
    avg_v3_eff = mean(v3_eff_0.22, na.rm = TRUE),
    .groups = "drop"  # Avoids keeping grouped structure in the result
  ) 
  average_efficiency$Location_Station <- paste(average_efficiency$Location, average_efficiency$Station_Number, sep = "_")




ggplot(average_efficiency, aes(x = Timepoint, y = avg_bac_eff)) +
  # geom_line(aes(group = interaction(Location_Station, Replicate)), size = 1) +
  geom_smooth(method = 'lm', se= F) +
  geom_point(size = 3) +
  facet_wrap(~ Location_Station, scales = "free") +
  labs(
    title = "0.22 µm Bacterial efficiency  Over Time",
    x = "Timepoint",
    y = "Viral Counts (c_Viruses)",
    color = "Replicate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )
