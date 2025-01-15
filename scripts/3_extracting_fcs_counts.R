
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
project_title<- "PE_Cruises_FCS"


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

metadata_processing(paste0(proj_dir,"/metadata/PE_Cruises_VP_Metadata.xlsx"), extension = ".xlsx", project_title = "PE_Cruises")

# #I will susbet metadata for testing rn
# metadata<- metadata[1:10,]

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


#As we see that gates don't fit all samples the same, we'll have to change the gating coordinates.
#You can use populate_gate_df function to do so

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
