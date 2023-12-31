#1.0 Function to install and load packages ####
install_and_load_packages <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("flowCore", "flowWorkspace", "ggcyto")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

# List of packages you want to install and load
packages_to_load <- c("tidyverse", 
                      "flowWorkspace",
                      "flowCore",
                      "scales",
                      "readxl",
                      "ggcyto",
                      "ggsci",
                      "svglite")

# Call the function with the list of packages
install_and_load_packages(packages_to_load)




#2.0 Setting up directories for data, metadata, results, and reference files####
set_up_vp_count<- function(){
  

  { #create directories
    dir_to_create<- c("/data", "/results", "/data/raw_data", 
                      "/data/metadata", "/data/ref_fcs")
    
    dir_to_create <- paste0(work_dir, dir_to_create)
    
    
    # Create the directories if they don't already exist
    for (dir in dir_to_create) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      } else {
        cat("Directory", dir, "already exists.\n")
      }
    }
    .GlobalEnv$work_dir <- paste0(getwd(), "/", project_title)
   
    rm(dir_to_create)    
        
  }
}



#3.0 Metadata processing
metadata_processing<- function(file_path, extension = ".xlsx", sheet =NULL, project_title, format = "%y%m%d", ...){
  if (extension == ".xlsx"){
    metadata<- read_excel(file_path, sheet)
  } else if (extension == ".csv") {
    metadata<- read.csv(file_path)
  } else {
    print("File format not supported. Please convert the metadata file to .xlsx or .csv")
  }
  
  metadata<- metadata[, c('Sample_Name', 'Staining_Protocol', 'Date_Measurement', 
                          'Location', 'Expt_Date', 'Station_Number', 'Depth', 
                          'Sample_Type', 'Timepoint', 'Replicate', 'Acquisition_Duration',
                          'Dilution', 'Flowrate')] %>%
    mutate(Sample_Index = row_number(),
           Expt_Date = as.Date(as.character(Expt_Date), format)) %>% #converting Expt_Date and Date_Measurement into Date format
    mutate(Date_Measurement= as.Date(as.character(Date_Measurement), format))
  .GlobalEnv$metadata<- metadata
  
  write.csv(metadata, file = paste0(work_dir, "./data/metadata/",project_title,"_metadata.csv"), row.names=F)
  print("Metadata processed and stored under `data/metadata`")
  return(metadata)
  
  
}


##Import doesn't work

import_fcs<- function(fcs_dir, project_title = "Project", ...){
  
  
  print("It'll take some time before all the files are transferred")
  for (sample in  1:length(metadata$Sample_Name)){
    file_source<- paste0(fcs_dir, metadata$Sample_Name[sample])
    file_dest <-  paste0(work_dir, "/data/raw_data")
    file.copy(from = file_source,
            to = file_dest,
            recursive = F)
  }
  if( "FALSE" %in% (metadata$Sample_Name %in% list.files("/data/raw_data/"))) { #files are missing in the raw_data folder
    print("Missing files. Check `missing_fcs_file`")
    print(setdiff(metadata$Sample_Name,list.files("/data/raw_data/")))#what to do with the missing files
    missing_fcs_file<- setdiff(metadata$Sample_Name,list.files("/data/raw_data/"))
    .GlobalEnv$missing_fcs_file <- missing_fcs_file
  } else {
    print("All files were transferred.")
  }
  fcs_data<- paste0(work_dir,"/data/raw_data/", 
                    list.files(paste0(work_dir, "/data/raw_data/")))
}




populate_gate_df <- function(sample_range = c(1: length(metadata$Sample_Index)),
                             g2_b_ssc = b_ssc, g2_b_fl1 = b_fl1, 
                             g2_v_ssc = v_ssc, g2_v_fl1 = v_fl1, 
                             g2_hna_ssc = hna_ssc, g2_hna_fl1 = hna_fl1, 
                             g2_lna_ssc = lna_ssc, g2_lna_fl1 = lna_fl1, 
                             g2_v1_ssc = v1_ssc, g2_v1_fl1 = v1_fl1, 
                             g2_v2_ssc = v2_ssc, g2_v2_fl1 = v2_fl1, 
                             g2_v3_ssc = v3_ssc, g2_v3_fl1 = v3_fl1,
                             adj_x_all_by = 0, adj_y_all_by = 0, VP_VPC_different = F) {
  
  
  
  if (VP_VPC_different == F){
    {
      b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6)
      b_fl1 = c(1.9, 1.4, 1.4,1.8, 2.8, 3.7, 3.2, 2.75)
      v_ssc = c(-0.25, 1.2)
      v_fl1 = c(-0.1, 1.7)
      hna_ssc = c(0.6, 3.5)
      hna_fl1 = c(2.15, 3.5)
      lna_ssc = c(1.0, 3.5)
      lna_fl1 = c(1.4, 2.15)  
      v1_ssc = c(-0.1, 0.90)
      v1_fl1 = c(-0.1, 0.8)
      v2_ssc = c(-0.1, 0.90)
      v2_fl1 = c(0.8, 1.25)
      v3_ssc = c(-0.1, 1.3)
      v3_fl1 = c(1.25, 1.7)
    }
  }else{
    {
      b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6)
      b_fl1 = c(1.9, 1.4, 1.4,1.8, 2.8, 3.7, 3.2, 2.75)
      v_ssc = c(-0.25, 1.2)
      v_fl1 = c(-0.1, 1.7)
      hna_ssc = c(0.6, 3.5)
      hna_fl1 = c(2.15, 3.5)
      lna_ssc = c(1.0, 3.5)
      lna_fl1 = c(1.4, 2.15)  
      v1_ssc = c(-0.1, 0.90)
      v1_fl1 = c(-0.1, 0.9)
      v2_ssc = c(-0.1, 0.90)
      v2_fl1 = c(0.9, 1.25)
      v3_ssc = c(-0.1, 1.3)
      v3_fl1 = c(1.25, 1.7)
    }
  }
  
  
  if (!exists("gate_df")) {
  gate_df<- data.frame(matrix(ncol = 16, nrow = nrow(metadata)))
  gate_df[,1:2]<- metadata[,c('Sample_Index', 'Sample_Name')]
  
  colname<- c("Sample_Index", "Sample_Name",
              "g_b_ssc", "g_b_fl1",
              "g_v_ssc", "g_v_fl1",
              "g_hna_ssc", "g_hna_fl1",
              "g_lna_ssc", "g_lna_fl1",
              "g_v1_ssc", "g_v1_fl1",
              "g_v2_ssc", "g_v2_fl1",
              "g_v3_ssc", "g_v3_fl1")
  colnames(gate_df)<- colname
  
  
  for(col in colname[-c(1:2)]) {
    gate_df[[col]] <- vector("list", nrow(gate_df))
  }
  }
  
  # Range of samples to update
  rows_to_update <- which(gate_df$Sample_Index %in% sample_range)
  
  adjssc <- 0
  adjfl1 <- 0 
  
  adj_args <- list(adj_x_all_by = adj_x_all_by, 
                   adj_y_all_by = adj_y_all_by)
  
  
  for (var_name in names(adj_args)) {
    cat("Variable", var_name, "has value:", adj_args[[var_name]], "\n")
    
    if (var_name == "adj_x_all_by") {
      cat("All SSC coordinates were adjusted by", adj_args[[var_name]], "\n")
      adjssc <- adjssc + adj_args[[var_name]]
    } 
    
    if (var_name == "adj_y_all_by") {
      cat("All FL1 coordinates were adjusted by", adj_args[[var_name]], "\n")
      adjfl1 <- adjfl1 + adj_args[[var_name]]
    } 
    
  }
  
  gate_df[rows_to_update, ] <- gate_df[rows_to_update, ] %>%
    mutate(
      g_b_ssc = list(list(g2_b_ssc )),
      g_b_fl1 = list(list(g2_b_fl1 )),
      g_v_ssc = list(list(g2_v_ssc )),
      g_v_fl1 = list(list(g2_v_fl1 )),
      g_hna_ssc = list(list(g2_hna_ssc + adjssc)),
      g_hna_fl1 = list(list(g2_hna_fl1 + adjfl1)),
      g_lna_ssc = list(list(g2_lna_ssc + adjssc)),
      g_lna_fl1 = list(list(g2_lna_fl1 + adjfl1)),
      g_v1_ssc = list(list(g2_v1_ssc + adjssc)),
      g_v1_fl1 = list(list(g2_v1_fl1 + adjfl1)),
      g_v2_ssc = list(list(g2_v2_ssc + adjssc)),
      g_v2_fl1 = list(list(g2_v2_fl1 + adjfl1)),
      g_v3_ssc = list(list(g2_v3_ssc + adjssc)),
      g_v3_fl1 = list(list(g2_v3_fl1 + adjfl1))
    )
  
  gate_df<<- gate_df
  return(gate_df)
}

gates<- function(gate_index, gate_df2 = gate_df, ... #index
                 ){
  polycut<- matrix(c(gate_df2$g_b_ssc[gate_index][[1]][[1]], gate_df2$g_b_fl1[gate_index][[1]][[1]]), 
                   nrow = length(gate_df2$g_b_ssc[gate_index][[1]][[1]]), 
                   ncol=2)
  colnames(polycut)<- c("SSC-H", "FL1-H")
  bgate<<- polygonGate(.gate = polycut, filterId = "Bacteria" )
  vgate<<- rectangleGate(filterId= "Viruses", 
                         "SSC-H" = gate_df2$g_v_ssc[gate_index][[1]][[1]], 
                         "FL1-H" = gate_df2$g_v_fl1[gate_index][[1]][[1]]) 
  
  HNA_Bacteria_bv<<- rectangleGate(filterId="HNA_Bacteria", 
                                   "SSC-H" = gate_df2$g_hna_ssc[gate_index][[1]][[1]], 
                                   "FL1-H" = gate_df2$g_hna_fl1[gate_index][[1]][[1]]) 
  LNA_Bacteria_bv<<- rectangleGate(filterId="LNA_Bacteria", 
                                   "SSC-H" = gate_df2$g_lna_ssc[gate_index][[1]][[1]], 
                                   "FL1-H" = gate_df2$g_lna_fl1[gate_index][[1]][[1]]) 
  
  #same viral gates as we don't utilise the viral info from bacterial samples
  v1<<- rectangleGate(filterId="V1", 
                      "SSC-H" = gate_df2$g_v1_ssc[gate_index][[1]][[1]], 
                      "FL1-H" = gate_df2$g_v1_fl1[gate_index][[1]][[1]])  
  v2<<- rectangleGate(filterId="V2", 
                      "SSC-H" = gate_df2$g_v2_ssc[gate_index][[1]][[1]], 
                      "FL1-H" = gate_df2$g_v2_fl1[gate_index][[1]][[1]]) 
  v3<<- rectangleGate(filterId="V3", 
                      "SSC-H" = gate_df2$g_v3_ssc[gate_index][[1]][[1]], 
                      "FL1-H" = gate_df2$g_v3_fl1[gate_index][[1]][[1]]) 
  
  detectors<<- c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H")
  
  translist_bv<<- transformList(detectors, logTransform())
  
}
  


ref_fcs_create<- function(ref_fcs_file){
  
  
  file.copy(from = paste0(work_dir,"/data/raw_data/", ref_fcs_file),
            to = paste0(work_dir,"/data/ref_fcs"))
  file.rename(from = paste0(work_dir,"/data/ref_fcs/", ref_fcs_file),
              to = paste0(work_dir,"/data/ref_fcs/", ref_fcs_file,"_ref"))
  ref_fcs<- paste0(work_dir, "/data/ref_fcs/", ref_fcs_file, "_ref")
  
  .GlobalEnv$ref_fcs<- ref_fcs
  print("Reference FCS file was created and is stored as:")
  print(ref_fcs)
  return(ref_fcs)
}

#Gating Functions
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

bacterial_gate<- function(gate = "same", df = metadata, fcs_file, ...){
  
  if (gate == "different"){
    if(df[df$Sample_Name == fcs_file,]$Staining_Protocol == "Bacteria") { #MAKE THIS INTO A FUNCTION
      print("Bacteria")
      .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_b
      .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_b
    } else if (df[df$Sample_Name == fcs_file,]$Staining_Protocol == "Viruses"){
      print("Viruses")
      .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_v
      .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_v
      
    } else {
      print("Staining protocol niether Bacteria or Viruses")
    }
  } else if (gate == "same"){
    .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_bv
    .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_bv
  }
  
  
}

read_transform_fs_bv <- function(x, ...){ #function to read and transform fcs files
  flowCore::read.flowSet(c(x, x)) %>%
    flowCore::transform(translist_bv)
  
}

gatingset_bv_stats<- function(gate_index, flowset, ...){ #flowset here is already transformed, cleaned, and compensated
  
  gates(gate_index)
  
  gsbv_fs<- flowWorkspace::GatingSet(flowset)
  gs_pop_add(gsbv_fs, bgate, parent="root")
  gs_pop_add(gsbv_fs, HNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, LNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, vgate, parent= "root")
  gs_pop_add(gsbv_fs, v1, parent="Viruses")
  gs_pop_add(gsbv_fs, v2, parent="Viruses")
  gs_pop_add(gsbv_fs, v3, parent="Viruses")
  recompute(gsbv_fs)
  stats<- gs_pop_get_stats(gsbv_fs[[2]], c("root", "Bacteria", "HNA_Bacteria",
                                           "LNA_Bacteria", "Viruses",
                                           "V1", "V2", "V3"), "count")
  write.table(stats, file = csv_file, sep = ",", append = T, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  
}


get_bv_stats<- function(df = metadata, gate = "same", write_csv = T, test = F, ...){
  
  if(test == T){
    df<- df[1:10,]
  }else if(test ==F){
    df<- df
  }
  
  
  if (write_csv == T){
    #create a csv file to store all the plots
    x<-data.frame(file_name = "file_name", pop = "pop", count = "count")
    x<-  x[-c(1),]
    
    if ("project_title" %in% ls(.GlobalEnv) == T){
      if(test == T){
        .GlobalEnv$project_title2<- paste0(project_title, "_test")
      }else if(test == F){ 
        .GlobalEnv$project_title2<- project_title
      }
      write.table(x,file= paste0(work_dir, "/results/",project_title2, "_counts.csv"), 
                  sep = ",", col.names = c("file_name", "pop", "count"))
      .GlobalEnv$csv_file<- paste0(work_dir, "/results/",project_title2, "_counts.csv")
      print(paste0(work_dir, "/results/",project_title,"_stats.csv created"))
    } else{
      write.table(x,file= paste0(work_dir, "/results/counts.csv"), 
                  sep = ",", col.names = c("file_name", "pop", "count"))
      .GlobalEnv$csv_file<- paste0(work_dir, "/results/counts.csv")
      print(paste0(work_dir, "/results/counts.csv", "created"))
    }
  }
  
  t<-0 #counter to 0
  
  
  
  for( i in 1:length(df$Sample_Index)){
    .GlobalEnv$i<- i
    .GlobalEnv$fcs_file<- df[df$Sample_Index == i,]$Sample_Name
    
    bacterial_gate(fcs_file = fcs_file)
    
    print(i)
    print(df[df$Sample_Index == i,]$Sample_Name)
    
    .GlobalEnv$fcs_data<- paste0(paste0(work_dir,"data/raw_data/", fcs_file))
    
    gates(gate_index = i)
    
    
    
    try(read_transform_fs_bv(fcs_data)  %>%
          gatingset_bv_stats(gate_index = i))
    t<- t +1
    print(paste0(t,"/",length(df$Sample_Name)))
    
    
    #print(HNA_Bacteria)
    #print(LNA_Bacteria)
    #rm(HNA_Bacteria, LNA_Bacteria, i, fcs_file, fcs_data)
  }
  
  # After the CSV has been written and before the function ends:
  if (write_csv == T){
    # Read the CSV back into R
    modified_df <- read.csv(.GlobalEnv$csv_file, stringsAsFactors = FALSE)
    
    # Remove .1 from Sample_Name values
    modified_df$file_name <- gsub("\\.1$", "", modified_df$file_name)
    
    # Write the modified dataframe back to the CSV
    write.csv(modified_df, .GlobalEnv$csv_file, row.names = FALSE)
  }
  
  rm(t)
  
  
}

gatingset_bv_plots<- function(gate_index, flowset, bins = 600, ...){ #flowset here is already transformed, cleaned, and compensated
  
  gates(gate_index)
  
  gsbv_fs<- flowWorkspace::GatingSet(flowset)
  gs_pop_add(gsbv_fs, bgate, parent="root")
  gs_pop_add(gsbv_fs, HNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, LNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, vgate, parent= "root")
  gs_pop_add(gsbv_fs, v1, parent="Viruses")
  gs_pop_add(gsbv_fs, v2, parent="Viruses")
  gs_pop_add(gsbv_fs, v3, parent="Viruses")
  recompute(gsbv_fs)
  p<- ggcyto::ggcyto(gsbv_fs[[2]], aes(x = `SSC-H`, y = `FL1-H`), subset = "root") +
    geom_hex(bins = bins, ... ) +  
    theme_bw()+
    labs(title= paste0(gate_index, "   ", metadata[metadata$Sample_Index== gate_index, ]$Sample_Name,  ...), x = "Side scatter (a.u.)", y = "Green Fluorescence (a.u.)")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "white")) +
    ggcyto_par_set(limits = list(x = c(-0.75,3.8), y = c(0, 3.6)))
  p<- p + geom_gate("Viruses", colour = "black", size = 1) + geom_stats("Viruses", type = c("gate_name", "count"), location = 'plot', adjust = c(0.9, 0.2), colour = "black")
  p<- p + geom_gate("Bacteria", colour = "black", size = 1) + geom_stats("Bacteria", type = c("gate_name", "count"), location = 'plot',  adjust = c(0.4, 1), colour = "black")
  p<- p + geom_gate(c("HNA_Bacteria", "LNA_Bacteria"), colour = "red", size = 0.8) + geom_stats(c("HNA_Bacteria", "LNA_Bacteria"), type = c("gate_name", "count"#, "percent"
                                                                                                                                            ), adjust = c(-0.4, 0.8), colour = "red")
  p<- p + geom_gate(c("V1", "V2", "V3"), colour = "blue", size = 0.8) + geom_stats(c("V1", "V2", "V3"), type = c("gate_name", "count"#, "percent"
                                                                                                                 ),  adjust = c(-0.5,0.7), colour = "blue")
  
  print(p)
  
}

get_bv_plots<- function(df = metadata, gate = "same", write_pdf = T, test = F, ...){
  
  if(test == T){
    df<- df[1:10,]
  }else if(test ==F){
    df<- df
  }
  
  
  #create a PDF file to store all the plots
  
    pdf(paste0(work_dir,"/results/",project_title,"_plots.pdf"), onefile =T)
    print(paste0(work_dir, "/results/",project_title,"_plots.pdf created"))
  
     
  if (write_pdf == T){
    
    
      if(test == T){
        .GlobalEnv$project_title2<- paste0(project_title, "_test")
      }else if(test == F){ 
        .GlobalEnv$project_title2<- project_title
      }
      pdf(paste0(work_dir,"/results/",project_title2,"_plots.pdf"), onefile =T)
      print(paste0(work_dir,"/results/",project_title2,"_plots.pdf created"))
      
     
  }
  
  t<-0 #counter to 0
  p<- list() #creating an empty list
  
  for( i in 1:length(df$Sample_Index)){
    .GlobalEnv$i<- i
    .GlobalEnv$fcs_file<- df[df$Sample_Index == i,]$Sample_Name
    
    bacterial_gate(fcs_file = fcs_file)
    
    print(i)
    print(df[df$Sample_Index == i,]$Sample_Name)
    
    .GlobalEnv$fcs_data<- paste0(paste0(work_dir,"/data/raw_data/", fcs_file))
    
    gates(i)
    
    p[[i]]<- list()
    p[[i]]<- try(read_transform_fs_bv(fcs_data)  %>%
                   gatingset_bv_plots(gate_index = i))
    t<- t +1
    print(paste0(t,"/",length(df$Sample_Index)))
    
    
    print(HNA_Bacteria)
    print(LNA_Bacteria)
    #rm(HNA_Bacteria, LNA_Bacteria, i, fcs_file, fcs_data)
  }
  rm(t)
  graphics.off()
  
  
}

cytoplot<- function(index, metadata_df = metadata, bins = 600, ...){
  
  fcs_file<- metadata_df[metadata_df$Sample_Index[index], ]$Sample_Name
  
  .GlobalEnv$fcs_data<- paste0(paste0(work_dir,"/data/raw_data/", fcs_file))
  
  plot<- read_transform_fs_bv(fcs_data, ...)  %>%
    gatingset_bv_plots(gate_index = index, bins, ...)
  
  return(plot)
  
}



combine_metadata_counts <- function(metadata_df = metadata, counts_df = counts){
  
  counts<- pivot_wider(counts_df, names_from = "pop", values_from = "count")
  #dim(counts)
  colnames(counts)[colnames(counts) == 'root'] <- "Total"
  colnames(counts)[colnames(counts) == 'file_name'] <- "Sample_Name"
  #head(counts)
  
  
  if (length(setdiff(counts$Sample_Name, metadata_df$Sample_Name)) > 0) {
    print("There are samples that were counted but doesn't have associated metadata")
  } else {
    print("All samples processed have associated metadata")
  }
  
  {
    counts_metadata<- base::merge(counts, metadata_df, by = "Sample_Name") 
    #counts_metadata<- dplyr::mutate(counts_metadata,Expt_Date = as.Date(as.character(Expt_Date), format= "%y%m%d"))
    #counts_metadata<- dplyr::mutate(counts_metadata,Date_Measurement = as.Date(as.character(Date_Measurement), format= "%y%m%d"))
    #flowrate<- dplyr::mutate(flowrate, Date_Measurement = as.Date(as.character(Date_Measurement), format= "%y%m%d"))
    #we could have done flowrate per date of measurement, but because one of our sample is from 2020, the flowrate wasn't noted. So let's go with average.
    # dim(counts_metadata)
    # head(counts_metadata)
    #counts_metadata[,'Flowrate'] <- mean(flowrate$Flow_rate)
  }
  
  #Adding Boolean columns for bacterial and viral total counts
  {
    counts_metadata$HNALNA <- counts_metadata$HNA_Bacteria + counts_metadata$LNA_Bacteria
    counts_metadata$V1V2V3 <- counts_metadata$V1 + counts_metadata$V2 + counts_metadata$V3
    
    #We can perform a simple linear regression to see if we could use Boolean 
    # addition of HNA/LNA and V1/V2/V3 as replacements for total bacterial and viral 
    # counts. Alternatively, these could have been added as Boolean gates during processing
    # scatter.smooth(counts_metadata$Bacteria, counts_metadata$HNALNA, main = "Bacterial Counts")
    #  summary(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA))  #Gives an R-square value of 1.0
    #Perhaps the next few lines aren't as important. 
    # plot(fitted(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)), resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    # abline(0,0)
    # qqnorm(resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    # qqline(resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    
    # scatter.smooth(counts_metadata$Viruses, counts_metadata$V1V2V3, main = "Viral Counts")
    #  summary(lm(counts_metadata$Viruses ~ counts_metadata$V1V2V3))  #Gives an R-square value of 0.999
    
  }
  .GlobalEnv$counts_metadata<- counts_metadata
  return(counts_metadata)
  
}

calc_TE<- function(df = counts_metadata){ 
  
  counts_metadata2<- df
  #
  #WHY?
  
  
  counts_metadata2[,'TE_Ba']<- NA
  counts_metadata2[,'TE_Vi']<- NA
  
  #We can now run TE calculation script
  for (name in counts_metadata2$Sample_Name){ #Here I added TE value for both viruses and bacteria
    if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Staining_Protocol == 'Viruses'
        #selecting rows with the specified file name that also had viral staining protocol
    ) {
      if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type  == 'TE') {
        print("TE") #we don't want this. so we move on and look for a count file
      } else if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type  != 'TE') {
        if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'
            #looking to see if there is a TE above this count file. if yes, we record it as 'a'
        ){
          a<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$HNALNA #HNALNA and not V1V2V3 because we might use viral samples for bacterial measurements in VP assay
          if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'
              #here we see if there is a TE above our first TE. we record this as 'b'
          ){
            b<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$HNALNA
          }
        }
        print(paste(name, a, b, mean(c(a,b)))) #too see what the output is
        counts_metadata2[counts_metadata2$Sample_Name == name,]$TE_Ba <- mean(c(a,b)) ##We are recording TE for bacteria in viral samples
        #adding the output to TE_value column
      }
      
      if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type  == 'TE') {
        print("TE") #we don't want this. so we move on and look for a count file
      } else if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type  != 'TE') {
        if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'
            #looking to see if there is a TE above this count file. if yes, we record it as 'a'
        ){
          c<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$V1V2V3
          if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'
              #here we see if there is a TE above our first TE. we record this as 'b'
          ){
            d<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$V1V2V3
          }
        }
        print(paste(name, a, b, mean(c(a,b)))) #too see what the output is
        counts_metadata2[counts_metadata2$Sample_Name == name,]$TE_Vi <- mean(c(c,d))
      }
    } else { #THIS PART WILL NOT RUN AS ALL MY SAMPLES ARE VIRAL 
      if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Staining_Protocol == 'Bacteria') {
        if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type  == 'TE') {
          print("yes")
        } else if (counts_metadata2[counts_metadata2$Sample_Name == name,]$Sample_Type != 'TE') {
          if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'){
            a<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-1), ]$HNALNA
            if (counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'){
              b<- counts_metadata2[which(counts_metadata2$Sample_Name == name) + c(-2), ]$HNALNA
            }
          }
          print(paste(name, a, b, mean(c(a,b))))
          counts_metadata2[counts_metadata2$Sample_Name == name,]$TE_Ba <- mean(c(a,b)) ##We are recording TE for bacteria in bacterial samples. There is no need to do it for viruses
        }
      }
    }
  }
  .GlobalEnv$counts_metadata_TE<- counts_metadata2
  return(counts_metadata2)
}

adjust_TE<- function(counts_metadata_df = counts_metadata, write_csv = T){
  
  counts_metadata2<- counts_metadata_df
  for (cols in c( "c_Bacteria", "c_HNA", "c_LNA", "c_Viruses", "c_V1", "c_V2", "c_V3")){
    counts_metadata2[, cols] <- NA
  }
  
  {
    counts_metadata2$c_Viruses<- with(counts_metadata2, ((V1V2V3-((V1V2V3/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$c_V1<- with(counts_metadata2, ((V1-((V1/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$c_V2<- with(counts_metadata2, ((V2-((V2/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$c_V3<- with(counts_metadata2, ((V3-((V3/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$c_Bacteria<- with(counts_metadata2, ((HNALNA-((HNALNA/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$c_HNA<- with(counts_metadata2, ((HNA_Bacteria-((HNA_Bacteria/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration)) 
    counts_metadata2$c_LNA<- with(counts_metadata2, ((LNA_Bacteria-((LNA_Bacteria/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata2$VBR<- with(counts_metadata2, c_Viruses/c_Bacteria) #calculated with bacteria from viral samples. Did an LR and the R2 was 0.95
    counts_metadata2$HNAperLNA<- with(counts_metadata2, c_HNA/c_LNA)
    #scatter.smooth(S22[S22$Staining_Protocol == 'Bacteria',]$c_Bacteria ~ S22[S22$Staining_Protocol == 'Viruses',]$c_Bacteria )
    #summary(lm(S22[S22$Staining_Protocol == 'Bacteria',]$c_Bacteria ~ S22[S22$Staining_Protocol == 'Viruses',]$c_Bacteria ))
  }
  
  {
    
    otpt_df<- counts_metadata2[counts_metadata2$Sample_Type != 'TE',]
    otpt_df<- otpt_df[, c('Sample_Name', 'Staining_Protocol', 'Expt_Date', 'Date_Measurement',
                          'Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', 'c_Bacteria', 'c_HNA', 'c_LNA', 
                          'c_Viruses', 'c_V1', 'c_V2', 'c_V3', 'VBR', 'HNAperLNA')]
    otpt_df<- otpt_df[
      with(otpt_df,
           order(
             otpt_df[, 'Staining_Protocol'],
             otpt_df[, 'Location'],
             otpt_df[, 'Station_Number'],
             otpt_df[, 'Depth'],
             otpt_df[, 'Timepoint'],
             otpt_df[, 'Replicate'],
             otpt_df[, 'Sample_Type']
             
           )),
    ]
  }
  if (write_csv == T){
    write.csv(otpt_df, paste0("PE_Cruises/results/", project_title, "_per_mL.csv"), row.names = F)
  }
  print("Counts were adjusted with TE")
  return(otpt_df)
  
}

bacterial_count_overview_plots<- function(data= NJ2020){
  plot1<- ggplot(data, aes(x = Sample_Type, y = c_Bacteria, col = as.factor(Station_Number), group = interaction(Sample_Type,Station_Number)))+
    geom_violin(aes(group = Sample_Type, fill = Sample_Type), alpha = 0.15, color = NA)+
    labs(fill = "Sample Type")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Bacterial Counts per Sample Type")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot1)
  ggsave("./results/TotalBacteria_perSampleType.svg", width = 20, height = 20, units = "cm")
  
  
  plot2<- ggplot(data, aes(x = as.factor(Station_Number), y = c_Bacteria, col = as.factor(Sample_Type), group = interaction(Sample_Type,Station_Number)))+
    geom_violin(aes(group = as.factor(Station_Number), fill = as.factor(Station_Number)), alpha = 0.25, color = NA)+
    labs(fill = "Expt No")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Bacterial Counts per Expt No")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot2)
  ggsave("./results/TotalBacteria_perExptNo.svg", width = 30, height = 20, units = "cm")
  
}

viral_count_overview_plots<- function(data = NJ2020){
  plot1<- ggplot(data, aes(x = Sample_Type, y = c_Viruses, col = as.factor(Station_Number), group = interaction(Sample_Type,Station_Number)))+
    geom_violin(aes(group = Sample_Type, fill = Sample_Type), alpha = 0.15, color = NA)+
    labs(fill = "Sample Type")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Viral Counts per Sample Type")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot1)
  ggsave("./results/TotalViruses_perSampleType.svg", width = 20, height = 20, units = "cm")
  
  plot2<- ggplot(data, aes(x = as.factor(Station_Number), y = c_Viruses, col = as.factor(Sample_Type), group = interaction(Sample_Type,Station_Number)))+
    geom_violin(aes(group = as.factor(Station_Number), fill = as.factor(Station_Number)), alpha = 0.25, color = NA)+
    labs(fill = "Expt No")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Viral Counts per Expt No")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  
  print(plot2)
  ggsave("./results/TotalViruses_perExptNo.svg", width = 30, height = 20, units = "cm")
  
}    
    
    

  
  
  
  
