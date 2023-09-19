
#The aim of this script is to process `fcs` files and metadata to extract per mL values.

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

#Step 1: 

#Define the project title

project_title<- "PE_Cruises"


# Define the working directory `work_dir`.\
# This could be done by using `getwd()` function.\

work_dir<- paste0(getwd(), "/", project_title, "/")


# **2. Setting up environment**\
# This is to create a space for your data, metadata, raw data, and ref_fcs file, 
# along with a directory for output.\
# Also, to install and load all the required libraries.\
# Use the function `set_up_vp_count` for this.\

source("scripts/1_importing_fcs_metadata_source.R")

#Setting up the directories using a project name
set_up_vp_count()

# **3.  Import and process metadata file**\
# The function `metadata_processing` accepts file path to the metadata file. It 
# also takes care of the dates.
# Please ensure that you have the following parameters mentioned in your metadata file.\
# ADD THE PARAMETERS HERE\

metadata_processing(paste0(getwd(),"/Metadata/PE_Cruises_VP_Metadata.xlsx"), extension = ".xlsx", project_title = "PE_Cruises")


# **4. Importing FCS files**\
# We can utilize the `Sample_Name` variable in the `metadata` file to copy files for analyses.\
# We need to provide the origin of these FCS files. Files will be copied to `./data/raw_data`.\
# 
# source_file_name <- metadata$Sample_Name[1]
# source_file_path <- file.path("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/PE_Cruises/data/raw_data", source_file_name)
# 
# if (file.exists(source_file_path)) {
#   cat("Source file exists.\n")
# } else {
#   cat("Source file does not exist.\n")
# }


for (sample in 1:length(metadata$Sample_Name)){
  print(metadata$Sample_Name[sample])
  file.copy(from = paste0(fcs_dir2, metadata$Sample_Name[sample]),
            to = paste0(work_dir, "/data/raw_data/"),
            recursive = T)
}

origin<- "C:/Users/hisham.shaikh/OneDrive - UGent/PE477_PE486_Viral_Production/Data/"
missing_fcs_file <- character(0)  # Initialize an empty vector for missing files

for (sample in seq_along(metadata$Sample_Name)) {
  source_file <- file.path(fcs_dir, metadata$Sample_Name[sample])
  destination_file <- file.path(work_dir, "data/raw_data", metadata$Sample_Name[sample])
  
  file.copy(from = source_file, to = destination_file)
  
  if (!file.exists(destination_file)) {
    missing_fcs_file <- c(missing_fcs_file, metadata$Sample_Name[sample])
  }
}



if (length(missing_fcs_file) > 0) {
  print("Missing files. Check `missing_fcs_file`")
  print(missing_fcs_file)
} else {
  print("All files were transferred.")
}




import_matching_files <- function(source_dir, dest_dir, metadata) {
  # Ensure trailing slashes in the directory paths
  source_dir <- file.path(source_dir)
  dest_dir <- file.path(dest_dir, "data/raw_data")
  
 
  # Get a list of sample names from metadata
  sample_names <- metadata$Sample_Name
  
  # Loop through each sample name and copy matching files
  for (sample_name in sample_names) {
    source_file <- file.path(source_dir, sample_name)
    dest_file <- file.path(dest_dir, sample_name)
    
    if (file.exists(source_file)) {
      file.copy(from = source_file, to = dest_file, overwrite = TRUE)
      cat("Copied file:", source_file, "to", dest_file, "\n")
    } else {
      cat("Source file not found:", source_file, "\n")
    }
  }
}

working_directory <- "C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/"
source_directory <- "C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/Data/"
destination_directory <- "C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/PE_Cruises/data/raw_data/"




import_matching_files(source_dir = source_directory, dest_dir = destination_directory, metadata = metadata)


#ref_fcs_create("vi201130.062")


#### Gating ####

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


populate_gate_df() #This creates a default gating dataframe for all your samples.
#Now let's evaluate how it performs on our samples

#We always only start with plots

#get_bv_plots()

#Note changes that need to be made. You can visuaize a cytoplot using sample index
#and chnage gate_df using populate_df to make changes. 

#1
#Samples1 - 51 need no changes
#let's visualize one of them now

cytoplot(index = 46)

#2
#Samples 52 - 152 need to be lowered by 0.25

cytoplot(index = 62, bins = 300)

populate_gate_df(sample_range = c(52:152),
                 VP_VPC_different = T,
                 adj_y_all_by = -0.25 #from default
)

cytoplot(index = 62, bins = 100)


#I'm finding sme issuses separating VP and VPC V1 and V2 here, but total viruses works fine


#3
#Samples 153 - 166 need no changes

cytoplot(index = 160, bins = 200)

#4
#Samples 167 - 253 need to be lowered by 0.3

cytoplot(index = 250, bins = 200)


populate_gate_df(sample_range = c(167:253),
                 adj_y_all_by = -0.3 #from default
)

cytoplot(index = 250, bins = 200, )








#As we see that gates don't fit all samples the same, we'll have to change the gating coordinates.
#You can use populate_gate_df function to do so



cytoplot(index = 7)


##### i wil make these changes between V1 and V2 later. Fo rnow, total viruses are alright. 
#Let's get the stats out

get_bv_stats()


#Import csv
counts<- read.csv("./PE_Cruises/results/PE_Cruises_counts.csv")
str(counts)


combine_metadata_counts(metadata_df = metadata, counts_df = counts)
calc_TE(df = counts_metadata)
adjust_TE(counts_metadata_df = counts_metadata_TE)


counts_per_mL<- read.csv("./PE_Cruises/results/PE_Cruises_per_mL.csv")
