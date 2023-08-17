







project_title<- "PE_Cruises"

work_dir<- paste0(getwd(), "/")




source(here("scripts", "1_importing_fcs_metadata_source.R"))

#Setting up the directories using a project name
set_up_vp_count(project_title = "PE_Cruises")


metadata_processing(paste0(work_dir,"/Metadata/PE_Cruises_VP_Metadata.xlsx"), extension = ".xlsx", project_title = "PE_Cruises")



#Importing FCS will be a bit of trouble, as they're all over the place
origin<- paste0(getwd(),"/Data/")

import_fcs(origin, project_title = project_title)

ref_fcs_create("vi201130.062")

get_bv_plots()
get_bv_stats()
