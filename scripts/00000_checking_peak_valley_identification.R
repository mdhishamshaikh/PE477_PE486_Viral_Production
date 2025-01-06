library(dplyr)
library(ggplot2)

# Step 1: Group by Location, Station_Number, and Timepoint, then average replicates for each treatment
averaged_data <- filtered_data %>%
  group_by(Location, Station_Number, Timepoint, Sample_Type) %>%
  summarise(mean_c_Viruses = mean(c_Viruses, na.rm = TRUE),
            sem_c_Viruses = sd(c_Viruses, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

# Step 2: Separate VP and VPC treatments and calculate the difference
wide_data <- averaged_data %>%
  pivot_wider(names_from = Sample_Type, 
              values_from = c(mean_c_Viruses, sem_c_Viruses), 
              names_sep = "_") %>%
  mutate(across(starts_with("sem_c_Viruses"), ~ replace_na(.x, 0)))

# Step 3: Calculate the difference between VPC and VP treatments and their standard error
wide_data <- wide_data %>%
  mutate(
    difference_c_Viruses = mean_c_Viruses_VPC - mean_c_Viruses_VP,  # Difference between VPC and VP
    difference_sem = sqrt((sem_c_Viruses_VPC)^2 + (sem_c_Viruses_VP)^2)  # Simplified error calculation
  )


wide_data

wide_data_PE477_3 <- wide_data %>%
  dplyr::filter(Location == "PE477",
                Station_Number == 3) %>%
  arrange(Timepoint)  # Ensure timepoints are in order

wide_data_PE477_1 <- wide_data %>%
  dplyr::filter(Location == "PE477",
                Station_Number == 1) %>%
  arrange(Timepoint)  # Ensure timepoints are in order
wide_data_PE486_7 <- wide_data %>%
  dplyr::filter(Location == "PE486",
                Station_Number == 7) %>%
  arrange(Timepoint)  # Ensure timepoints are in order


DF2 <-wide_data_PE477_3
DF2 <-wide_data_PE477_1

DF2 <-wide_data_PE486_7
# Sample DataFrame to test the function
index_peaks <- vp_determine_peaks(c(+10e+10, DF2$difference_c_Viruses , -10e+10))
index_valleys <- vp_determine_valleys(c(+10e+10, DF2$difference_c_Viruses , -10e+10))
index_peaks
index_valleys
# Peak at 5 and valley at 3

index_peaks <- vp_determine_peaks_with_se(c(+10e+10, DF2$difference_c_Viruses, -10e+10),
                                          c(0, DF2$difference_sem , 0))
index_valleys <- vp_determine_valleys_with_se(c(+10e+10, DF2$difference_c_Viruses, -10e+10),
                                              c(0, DF2$difference_sem , 0))

index_peaks
index_valleys
# Fails for PE477_1. Couldn't recognize between T6 and T12 cause it only relies on neighbors

# Trying pracma
library(pracma)

# Wihtout SE
counts <- c(10e10, DF2$difference_c_Viruses, -10e10)

# Find peaks
findpeaks(counts)[,2]-1
#Find valleys
findpeaks(-counts)[,2]-1



#SEM 
sem <- c(0, DF2$difference_sem, 0)

peak_indices <- findpeaks(counts)
valley_indices <- findpeaks(neg_counts)

# Extract positions of peaks and valleys
peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys


valid_pairs <- list()
for (valley_pos in valley_positions) {
  # Find the first peak that comes after the valley
  next_peak_index <- which(peak_positions > valley_pos)[1]
  
  # Skip if no peak is found to the right
  if (is.na(next_peak_index)) next
  
  # Get the position of this next peak
  peak_pos <- peak_positions[next_peak_index]
  
  # Calculate means and SEMs
  valley_mean <- counts[valley_pos]
  valley_sem <- sem[valley_pos]
  
  peak_mean <- counts[peak_pos]
  peak_sem <- sem[peak_pos]
  
  # Step 3: Check SEM overlap condition
  if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
    # If no overlap, keep the pair
    valid_pairs[[length(valid_pairs) + 1]] <- c(valley_pos, peak_pos)
  }
}




# Load the necessary library
library(pracma)

# Your data: Assuming DF2$difference_c_Viruses and DF2$difference_sem are your numeric vectors
counts <- DF2$difference_c_Viruses
sem <- DF2$difference_sem

# Step 1: Add a large constant to shift all values to positive range
shift_constant <- abs(min(counts)) + 1e5  # Make sure all values are positive
shifted_counts <- counts + shift_constant

# Also shift the SEM to maintain the correct relative difference
shifted_counts_with_bounds <- c(shift_constant * 2, shifted_counts, shift_constant * -2)
shifted_sem <- c(0, sem, 0)

# Step 2: Identify peaks and valleys based on the shifted counts
peak_indices <- findpeaks(shifted_counts_with_bounds)
valley_indices <- findpeaks(-shifted_counts_with_bounds)  # Use negative to find valleys

# Extract positions of peaks and valleys
peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys

# Initialize a list to store valid valley-peak pairs
valid_pairs <- list()

# Step 3: Pair each valley with the nearest peak to the right
for (valley_pos in valley_positions) {
  # Find the first peak that comes after the valley
  next_peak_index <- which(peak_positions > valley_pos)[1]
  
  # Skip if no peak is found to the right
  if (is.na(next_peak_index)) next
  
  # Get the position of this next peak
  peak_pos <- peak_positions[next_peak_index]
  
  # Calculate means and SEMs for the original (non-shifted) data
  valley_mean <- counts[valley_pos - 1]  # Adjust position due to boundary added
  valley_sem <- sem[valley_pos - 1]
  
  peak_mean <- counts[peak_pos - 1]  # Adjust position due to boundary added
  peak_sem <- sem[peak_pos - 1]
  
  # Step 4: Check SEM overlap condition
  if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
    # If no overlap, keep the pair
    valid_pairs[[length(valid_pairs) + 1]] <- c(valley_pos - 1, peak_pos - 1)  # Correct for boundary
  }
}

# Output valid valley-peak pairs as a data frame
valid_pairs_df <- do.call(rbind, valid_pairs)
colnames(valid_pairs_df) <- c("Valley Position", "Peak Position")

# Print the confirmed valley-peak pairs
print(valid_pairs_df)



















###########################
counts <- c(10e10, DF2$difference_c_Viruses, -10e10) 
sem <- c(0, DF2$difference_sem, 0)




vp_determine_peaks_with_se <- function(counts, sem) {
  # Load the necessary library
  library(pracma)
  
  # Step 1: Add a large constant to shift all values to positive range
  shift_constant <- abs(min(counts)) + 10e10  # Ensure all values are positive
  shifted_counts <- counts + shift_constant
  
  # Step 2: Identify peaks and valleys in the shifted data
  peak_indices <- findpeaks(shifted_counts)
  valley_indices <- findpeaks(-shifted_counts)  # Use negative to find valleys
  
  # Extract positions of peaks and valleys
  peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
  valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys
  
  # Initialize vectors for valid indices
  valid_peaks <- c()
  valid_valleys <- c()
  
  # Step 3: Pair each valley with the nearest peak to the right and validate with SEM
  for (valley_pos in valley_positions) {
    # Find the first peak that comes after the valley
    next_peak_index <- which(peak_positions > valley_pos)[1]
    
    # Skip if no peak is found to the right
    if (is.na(next_peak_index)) next
    
    # Get the position of this next peak
    peak_pos <- peak_positions[next_peak_index]
    
    # Calculate means and SEMs for the original (non-shifted) data
    valley_mean <- shifted_counts[valley_pos]  # Adjust position due to boundary
    valley_sem <- sem[valley_pos]
    
    peak_mean <- shifted_counts[peak_pos]  # Adjust position due to boundary
    peak_sem <- sem[peak_pos]
    
    # Step 4: Check SEM overlap condition
    if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
      # If no overlap, store the valid indices
      valid_peaks <- c(valid_peaks, peak_pos - 1)
      valid_valleys <- c(valid_valleys, valley_pos - 1)
    }
  }
  
  # Return a list containing vectors of peak and valley indices
  return(list(peaks = valid_peaks, valleys = valid_valleys))
}


















determine_vp_VIPCAL_average_replicates_SE <- function(AVG_dataframe) {
  result_list <- list()
  
  for (combi_tag in unique(AVG_dataframe$tag)) {
    for (virus in unique(AVG_dataframe[AVG_dataframe$Microbe == 'Viruses', ]$Population)) {
      for (sample in unique(AVG_dataframe$Sample_Type)) {
        DF <- AVG_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Population == virus, .data$Sample_Type == sample)
        
        for (time in unique(DF$Time_Range)) {
          DF2 <- DF %>%
            dplyr::filter(.data$Time_Range == time)
          
          # Use the new function to get indices of peaks and valleys
          indices <- vp_determine_peaks_with_se(c(+10e+10, DF2$Mean, -10e+10),
                                                c(0, DF2$SE, 0))
          
          index_peaks <- indices$peaks
          index_valleys <- indices$valleys
          
          if (length(index_peaks) == 0) {
            viral_production <- 0
            abs_vp <- 0
            se <- 0
          } else {
            total_vp <- 0
            total_abs_vp <- 0
            total_se <- 0
            
            for (index in 1:length(index_peaks)) {
              viral_production_index <- (DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]) / 
                (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
              total_vp <- total_vp + viral_production_index
              
              abs_vp_index <- DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]
              total_abs_vp <- total_abs_vp + abs_vp_index
              
              se_index <- (DF2$SE[index_peaks[index]] + DF2$SE[index_valleys[index]]) / 
                (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
              total_se <- total_se + se_index
            }
            viral_production <- total_vp / length(index_peaks)
            abs_vp <- total_abs_vp
            se <- total_se / length(index_peaks)
          }
          result <- c(combi_tag, time, virus, sample, viral_production, abs_vp, se)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE')
  viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(AVG_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}








######################


library(dplyr)

# Define the function to determine peaks and valleys with SEM consideration
vp_determine_peaks_with_se <- function(counts, sem) {
  # Load necessary library
  library(pracma)
  
  # Step 1: Add a large constant to shift all values to positive range
  shift_constant <- abs(min(counts)) + 10e10  # Ensure all values are positive
  shifted_counts <- counts + shift_constant
  
  # Step 2: Identify peaks and valleys in the shifted data
  peak_indices <- findpeaks(shifted_counts)
  valley_indices <- findpeaks(-shifted_counts)  # Use negative to find valleys
  
  # Extract positions of peaks and valleys
  peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
  valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys
  
  # Initialize vectors for valid indices
  valid_peaks <- c()
  valid_valleys <- c()
  
  # Step 3: Pair each valley with the nearest peak to the right and validate with SEM
  for (valley_pos in valley_positions) {
    # Find the first peak that comes after the valley
    next_peak_index <- which(peak_positions > valley_pos)[1]
    
    # Skip if no peak is found to the right
    if (is.na(next_peak_index)) next
    
    # Get the position of this next peak
    peak_pos <- peak_positions[next_peak_index]
    
    # Calculate means and SEMs for the original (non-shifted) data
    valley_mean <- counts[valley_pos]  # Adjust position due to boundary
    valley_sem <- sem[valley_pos]
    
    peak_mean <- counts[peak_pos]  # Adjust position due to boundary
    peak_sem <- sem[peak_pos]
    
    # Step 4: Check SEM overlap condition
    if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
      # If no overlap, store the valid indices
      valid_peaks <- c(valid_peaks, peak_pos - 1)
      valid_valleys <- c(valid_valleys, valley_pos - 1)
    }
  }
  
  # Return a list containing vectors of peak and valley indices
  return(list(peaks = valid_peaks, valleys = valid_valleys))
}

# Main function to calculate viral production
determine_vp_VIPCAL_average_replicates_SE <- function(wide_data) {
  result_list <- list()
  
  # Loop over each unique Location and Station_Number
  for (location in unique(wide_data$Location)) {
    for (station in unique(wide_data$Station_Number)) {
      # Filter data for current Location and Station
      DF <- wide_data %>%
        filter(.data$Location == location, .data$Station_Number == station)
      
      # Use the new function to get indices of peaks and valleys
      indices <- vp_determine_peaks_with_se(c(+10e+10, DF$mean_c_Viruses_VP, -10e+10),
                                            c(0, DF$sem_c_Viruses_VP, 0))
      
      index_peaks <- indices$peaks
      index_valleys <- indices$valleys
      
      # Check if peaks and valleys are found
      if (length(index_peaks) == 0) {
        viral_production <- 0
        abs_vp <- 0
        se <- 0
      } else {
        total_vp <- 0
        total_abs_vp <- 0
        total_se <- 0
        
        for (index in 1:length(index_peaks)) {
          # Calculate viral production, absolute difference, and standard error
          viral_production_index <- (DF$mean_c_Viruses_VP[index_peaks[index]] - 
                                       DF$mean_c_Viruses_VP[index_valleys[index]]) / 
            (DF$Timepoint[index_peaks[index]] - DF$Timepoint[index_valleys[index]])
          total_vp <- total_vp + viral_production_index
          
          abs_vp_index <- DF$mean_c_Viruses_VP[index_peaks[index]] - DF$mean_c_Viruses_VP[index_valleys[index]]
          total_abs_vp <- total_abs_vp + abs_vp_index
          
          se_index <- (DF$sem_c_Viruses_VP[index_peaks[index]] + DF$sem_c_Viruses_VP[index_valleys[index]]) / 
            (DF$Timepoint[index_peaks[index]] - DF$Timepoint[index_valleys[index]])
          total_se <- total_se + se_index
        }
        
        # Average results
        viral_production <- total_vp / length(index_peaks)
        abs_vp <- total_abs_vp
        se <- total_se / length(index_peaks)
      }
      
      # Store the results
      result <- c(location, station, viral_production, abs_vp, se)
      result_list[[length(result_list) + 1]] <- result
    }
  }
  
  # Convert results to a data frame
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('Location', 'Station_Number', 'VP', 'abs_VP', 'VP_SE')
  viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')], as.numeric)
  
  return(viral_production_VIPCAL)
}

# Run the function on the wide_data
viral_production_results <- determine_vp_VIPCAL_average_replicates_SE(wide_data)

# Print the results
print(viral_production_results)




determine_vp_lysogenic__VIPCAL_average_replicates_SE <- function(wide_data) {
  result_list <- list()
  
  # Loop over each unique Location and Station_Number
  for (location in unique(wide_data$Location)) {
    for (station in unique(wide_data$Station_Number)) {
      # Filter data for current Location and Station
      DF <- wide_data %>%
        filter(.data$Location == location, .data$Station_Number == station)
      
      # Use the new function to get indices of peaks and valleys
      indices <- vp_determine_peaks_with_se(c(+10e+10, DF$difference_c_Viruses , -10e+10),
                                            c(0, DF$difference_sem , 0))
      
      index_peaks <- indices$peaks
      index_valleys <- indices$valleys
      
      # Check if peaks and valleys are found
      if (length(index_peaks) == 0) {
        viral_production <- 0
        abs_vp <- 0
        se <- 0
      } else {
        total_vp <- 0
        total_abs_vp <- 0
        total_se <- 0
        
        for (index in 1:length(index_peaks)) {
          # Calculate viral production, absolute difference, and standard error
          viral_production_index <- (DF$difference_c_Viruses [index_peaks[index]] - 
                                       DF$difference_c_Viruses [index_valleys[index]]) / 
            (DF$Timepoint[index_peaks[index]] - DF$Timepoint[index_valleys[index]])
          total_vp <- total_vp + viral_production_index
          
          abs_vp_index <- DF$difference_c_Viruses [index_peaks[index]] - DF$difference_c_Viruses [index_valleys[index]]
          total_abs_vp <- total_abs_vp + abs_vp_index
          
          se_index <- (DF$difference_sem [index_peaks[index]] + DF$difference_sem [index_valleys[index]]) / 
            (DF$Timepoint[index_peaks[index]] - DF$Timepoint[index_valleys[index]])
          total_se <- total_se + se_index
        }
        
        # Average results
        viral_production <- total_vp / length(index_peaks)
        abs_vp <- total_abs_vp
        se <- total_se / length(index_peaks)
      }
      
      # Store the results
      result <- c(location, station, viral_production, abs_vp, se)
      result_list[[length(result_list) + 1]] <- result
    }
  }
  
  # Convert results to a data frame
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('Location', 'Station_Number', 'VP', 'abs_VP', 'VP_SE')
  viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')], as.numeric)
  
  return(viral_production_VIPCAL)
}

# Run the function on the wide_data
viral_lysogenic_production_results <- determine_vp_lysogenic__VIPCAL_average_replicates_SE(wide_data)

# Print the results
print(viral_lysogenic_production_results)


viral_production_vp <- viral_production_results %>%
  left_join(viral_lysogenic_production_results, by = c("Location", "Station_Number"), suffix = c("_Lytic", "_Lysogenic"))

viral_production_vp$Depth <- 1

write.csv(viral_production_vp, file = "./results/viral_production_vipcal_se.csv", row.names = F)
