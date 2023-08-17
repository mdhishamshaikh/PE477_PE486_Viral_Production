
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

work_dir<- paste0(getwd(), "/../")


# **2. Setting up environment**\
# This is to create a space for your data, metadata, raw data, and ref_fcs file, 
# along with a directory for output.\
# Also, to install and load all the required libraries.\
# Use the function `set_up_vp_count` for this.\

source("scripts/1_importing_fcs_metadata_source.R")

#Setting up the directories using a project name
set_up_vp_count(project_title = project_title)

# **3.  Import and process metadata file**\
# The function `metadata_processing` accepts file path to the metadata file. It 
# also takes care of the dates.
# Please ensure that you have the following parameters mentioned in your metadata file.\
# ADD THE PARAMETERS HERE\

metadata_processing(paste0(work_dir,"/../Metadata/PE_Cruises_VP_Metadata.xlsx"), extension = ".xlsx", project_title = "PE_Cruises")





#Importing FCS will be a bit of trouble, as they're all over the place
origin<-"/../Data/"

import_fcs(origin, project_title = project_title)

ref_fcs_create("vi201130.062")

get_bv_plots()
get_bv_stats()
