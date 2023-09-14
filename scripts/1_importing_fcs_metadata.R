
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

source_file_name <- metadata$Sample_Name[1]
source_file_path <- file.path("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/Data/", source_file_name)

if (file.exists(source_file_path)) {
  cat("Source file exists.\n")
} else {
  cat("Source file does not exist.\n")
}


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

fcs_data <- list.files(file.path(work_dir, "data/raw_data"))

return(list(fcs_data = fcs_data, missing_fcs_file = missing_fcs_file))
}
import_fcs(origin)

wafile.exists("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/PE_Cruises/data/raw_data/vi201123.059")



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


ref_fcs_create("vi201130.062")
project_title<- "PE_Cruises_2"

#This is already in the source script, but if you want to adjust, you can do it here.
{
  polycut<- matrix(c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6,
                     1.9, 1.5, 1.5, 1.8, 2.8, 3.7, 3.2, 2.75), nrow = 8, ncol=2)
  colnames(polycut)<- c("SSC-H", "FL1-H")
  bgate<- polygonGate(.gate = polycut, filterId = "Bacteria" )
  vgate<- rectangleGate(filterId= "Viruses", "SSC-H" = c(-0.1,1.35), "FL1-H" = c(-0.1, 1.7)) 
  #Different bacterial gates for viral vs bacterial stained samples. We need different gates cause sometimes the bacterial
  #and viral samples have a slight shift.
  HNA_Bacteria_b<- rectangleGate(filterId="HNA_Bacteria",
                                 "SSC-H"=c(0.6, 3.5), "FL1-H"=c(2.15, 3.5))
  LNA_Bacteria_b<- rectangleGate(filterId="LNA_Bacteria",
                                 "SSC-H"=c(1.0, 3.5), "FL1-H"=c(1.5, 2.15))
  
  HNA_Bacteria_v<- rectangleGate(filterId="HNA_Bacteria",
                                 "SSC-H"=c(0.6, 3.5), "FL1-H"=c(2.3, 3.5))
  LNA_Bacteria_v<- rectangleGate(filterId="LNA_Bacteria",
                                 "SSC-H"=c(1.0, 3.5), "FL1-H"=c(1.5, 2.3))
  
  #in case it doesn't, we could define the gate as following
  HNA_Bacteria_bv<- rectangleGate(filterId="HNA_Bacteria",
                                  "SSC-H"=c(0.6, 3.7), "FL1-H"=c(2.3, 3.6))
  LNA_Bacteria_bv<- rectangleGate(filterId="LNA_Bacteria",
                                  "SSC-H"=c(1.0, 3.7), "FL1-H"=c(1.5, 2.3))
  
  #same viral gates as we don't utilise the viral info from bacterial samples
  v1<- rectangleGate(filterId="V1", 
                     "SSC-H"= c(-0.1, 0.90), "FL1-H"= c(-0.1, 0.8)) 
  v2<- rectangleGate(filterId="V2", 
                     "SSC-H"= c(-0.1, 0.90), "FL1-H"= c(0.8, 1.25)) 
  v3<- rectangleGate(filterId="V3", 
                     "SSC-H"= c(-0.1, 1.3), "FL1-H"= c(1.25, 1.7))
  
  detectors<- c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H")
  
  translist_bv<- transformList(detectors, logTransform())
}

gate_df<- data.frame(matrix(ncol = 29, nrow = nrow(metadata) ))
gate_df[,1]<- metadata[,1]

colnames(gate_df)<- c("Sample_Name", 
                      "v_x_s", "v_y_s",
                      "v_x_f", "v_y_f",
                      "b_x_s", "b_y_s",
                      "b_x_f", "b_y_f",
                      "v1_x_s", "v1_y_s",
                      "v1_x_f", "v1_y_f",
                      "v2_x_s", "v2_y_s",
                      "v2_x_f", "v2_y_f",
                      "v3_x_s", "v3_y_s",
                      "v3_x_f", "v3_y_f",
                      "hna_x_s", "hna_y_s",
                      "hna_x_f", "hna_y_f",
                      "lna_x_s", "lna_y_s",
                      "lna_x_f", "lna_y_f")

#provide an index for sample. Range
r1<- 1
r2<- 10

metadata = metadata
b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6)
b_fl1 = c(1.9, 1.5, 1.5, 1.8, 2.8, 3.7, 3.2, 2.75)
v_ssc = c(0.6, 3.5)
v_fl1 = c(-0.1, 1.7)
hna_ssc = c(0.6, 3.5)
hna_fl1 = c(2.15, 3.5)
lna_ssc = c(1.0, 3.5)
lna_fl1 = c(1.5, 2.15)
v1_ssc = c(-0.1, 0.90)
v1_fl1 = c(-0.1, 0.8)
v2_ssc = c(-0.1, 0.90)
v2_fl1 = c(0.8, 1.25)
v3_ssc = c(-0.1, 1.3)
v3_fl1 = c(1.25, 1.7)


  
  
gate_df[r1:r2, -c(1)]<- c(1,2,3)
a

coords<- c()

gate_df<- rbind(gate_df, c(as.character(metadata[i,1]), coords))



populate_gates<- function(){}


tib<- tibble(
  Sample_Name = as.character(metadata[i,1]),
  gate = list(
    tibble(b = v = tibble(x = v_x,
                      y = v_y),
           )
  )
)


get_bv_plots()

get_bv_stats()




##### plot for shiny

read_transform_fs_bv <- function(x){ #function to read and transform fcs files
  flowCore::read.flowSet(c(x, x)) %>%
    flowCore::transform(translist_bv)
  
}