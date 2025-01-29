# AIM: Extract viral decay rates

# 0.0 Setting up ####
source("scripts/0_source.R")

# 1.0 Importing FCS counts ####
counts_per_mL<- read.csv("./PE_Cruises_FCS/results/PE_Cruises_FCS_per_mL.csv")
str(counts_per_mL)

# Extracting 0.22 µm treatment
counts_0.22<- counts_per_mL %>%
  dplyr::filter(Sample_Type == '0.22',
                Staining_Protocol == 'Viruses') %>%
  mutate(Depth = ifelse(Depth == 1, 7, Depth)) %>%
  select(Location, Station_Number, Depth, Sample_Type, Timepoint, Replicate, c_Viruses, c_V1, c_V2, c_V3) %>%
  pivot_longer(cols = c(c_Viruses, c_V1, c_V2, c_V3),
               names_to = "Population",
               values_to = "Counts") %>%
  group_by(Location, Station_Number, Depth, Sample_Type, Timepoint, Population) %>%
  summarise(Mean_Counts = mean(Counts, na.rm = TRUE), .groups = 'drop')

# 2.0 Calculating decay rates ####

viral_decay <- counts_0.22 %>%
  pivot_wider(names_from = "Timepoint",
              values_from = "Mean_Counts") %>%
  mutate(decay_rate_linear = -(`24` - `0`)/24, # linear #making it positive
         decay_rate_exponential = (log(`0`) - log(`24`)) / 24) # exponential

# Adding original viral abundance to calculate percent viruses lost per day for linear
abundance_ts <- read.csv("./Data/nutrients_ts_abundance_coordinates/nutrients_ts_abundance_coordinates.csv")
abundance <- abundance_ts %>%
  select(Location, Station_Number, Depth, Total_Viruses, V1, V2, V3) %>%
  rename(c_Viruses = Total_Viruses,
         c_V1 = V1,
         c_V2 = V2,
         c_V3 = V3) %>%
  pivot_longer(cols = c(c_Viruses, c_V1, c_V2, c_V3),
               names_to = "Population",
               values_to = "Original_Abundance")

viral_decay_abundance<- viral_decay %>%
  left_join(abundance, by = c("Location", "Station_Number", "Depth", "Population"))

viral_decay_percent <- viral_decay_abundance %>%
  mutate(percent_decay_day_linear = (decay_rate_linear/Original_Abundance)*100*24, # percent per day for linear
         percent_decay_day_exponential = decay_rate_exponential *100*24) %>% # percent per day for exponential
  select(-`0`, -`24`, -Original_Abundance) 


write.csv(viral_decay_percent, "./results/viral_decay.csv", row.names = F, quote = F)


viral_decay_percent <- viral_decay_percent %>%
  mutate(combi_tag = paste(Location, Station_Number, sep = "_"))


viral_decay_long <- viral_decay_percent %>%
  tidyr::pivot_longer(
    cols = c(
      decay_rate_linear,
      decay_rate_exponential,
      percent_decay_day_linear,
      percent_decay_day_exponential
    ),
    names_to = c("Decay_Metric", "Decay_Type"),
    names_pattern = "(.*)_(linear|exponential)",
    values_to = "Value"
  ) %>%
  mutate(
    Decay_Metric = case_when(
      Decay_Metric == "decay_rate" ~ "Rate",
      Decay_Metric == "percent_decay_day" ~ "Percent",
      TRUE ~ Decay_Metric
    ),
    Decay_Type = case_when(
      Decay_Type == "linear" ~ "Linear",
      Decay_Type == "exponential" ~ "Exponential",
      TRUE ~ Decay_Type
    )
  )


# Plot
viral_decay_plot <- ggplot(viral_decay_long %>% dplyr::filter(Population == 'c_Viruses'), aes(x = combi_tag, y = Value, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = 'black') +
  facet_wrap(Decay_Type ~ Decay_Metric, scales = "free") +
  labs(x = "Stations", y = "Value", title = "Decay Rate and Percent Decay") +
  theme_test() +
  #scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
viral_decay_plot
ggsave(plot = viral_decay_plot, filename = "./figures/viral_decay_plot.svg", width = 12, height = 8, dpi = 400)



####################


# 1.0 Importing FCS counts ####
counts_per_mL<- read.csv("PE_Cruises_FCS_with_0.22/results/PE_Cruises_FCS_with_0.22_per_mL.csv")
# Creating a combined variable for Location and Station_Number
counts_per_mL$Location_Station <- paste(counts_per_mL$Location, counts_per_mL$Station_Number, sep = "_")


# 2.0 Visualizing 0.22 µm counts over time ####

# removing outliers
counts_per_mL_outliers_removed <- counts_per_mL %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE477" & Station_Number == 2 & Timepoint == 6 & Replicate == 1,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  )) %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE477" & Station_Number == 6 & Timepoint == 0 & Replicate == 2,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  )) %>%
  dplyr::mutate(c_Viruses = if_else(
    Location == "PE486" & Station_Number == 5 & Timepoint == 9 & Replicate == 3,
    NA_real_,  # Replace with NA for numeric columns
    c_Viruses
  ))


str(counts_per_mL_outliers_removed)

# Extracting 0.22 µm treatment
counts_0.22<- counts_per_mL_outliers_removed %>%
  dplyr::filter(Sample_Type == '0.22',
                Staining_Protocol == 'Viruses') %>%
  mutate(Depth = ifelse(Depth == 1, 7, Depth)) %>%
  select(Location, Station_Number, Depth, Sample_Type, Timepoint, Replicate, c_Viruses) %>%
  pivot_longer(cols = c(c_Viruses),
               names_to = "Population",
               values_to = "Counts") %>%
  group_by(Location, Station_Number, Depth, Sample_Type, Timepoint, Population) %>%
  summarise(Mean_Counts = mean(Counts, na.rm = TRUE), .groups = 'drop')

# 2.0 Calculating decay rates ####

viral_decay <- counts_0.22 %>%
  pivot_wider(names_from = "Timepoint",
              values_from = "Mean_Counts",
              names_prefix = "T" ) %>%
  mutate(decay_rate_linear = -(T24 - T0)/24, # linear #making it positive
         decay_rate_exponential = (log(T0) - log(T24)) / 24) # exponential

# Adding original viral abundance to calculate percent viruses lost per day for linear
abundance_ts <- read.csv("./Data/nutrients_ts_abundance_coordinates/nutrients_ts_abundance_coordinates.csv")
abundance <- abundance_ts %>%
  mutate(Depth = ifelse(Depth == 1, 7, Depth)) %>%
  select(Location, Station_Number, Depth, Total_Viruses, V1, V2, V3) %>%
  rename(c_Viruses = Total_Viruses,
         c_V1 = V1,
         c_V2 = V2,
         c_V3 = V3) %>%
  pivot_longer(cols = c(c_Viruses, c_V1, c_V2, c_V3),
               names_to = "Population",
               values_to = "Original_Abundance")

viral_decay_abundance<- viral_decay %>%
  left_join(abundance, by = c("Location", "Station_Number", "Depth", "Population"))

viral_decay_percent <- viral_decay_abundance %>%
  mutate(percent_decay_day_linear = (decay_rate_linear/Original_Abundance)*100*24, # percent per day for linear
         percent_decay_day_exponential = decay_rate_exponential *100*24) %>% # percent per day for exponential
  select(-c("T3", "T6", "T9", "T12" , "Original_Abundance")) 


write.csv(viral_decay_percent, "./results/viral_decay.csv", row.names = F, quote = F)


viral_decay_percent <- viral_decay_percent %>%
  mutate(combi_tag = paste(Location, Station_Number, sep = "_"))


viral_decay_long <- viral_decay_percent %>%
  tidyr::pivot_longer(
    cols = c(
      decay_rate_linear,
      decay_rate_exponential,
      percent_decay_day_linear,
      percent_decay_day_exponential
    ),
    names_to = c("Decay_Metric", "Decay_Type"),
    names_pattern = "(.*)_(linear|exponential)",
    values_to = "Value"
  ) %>%
  mutate(
    Decay_Metric = case_when(
      Decay_Metric == "decay_rate" ~ "Rate",
      Decay_Metric == "percent_decay_day" ~ "Percent",
      TRUE ~ Decay_Metric
    ),
    Decay_Type = case_when(
      Decay_Type == "linear" ~ "Linear",
      Decay_Type == "exponential" ~ "Exponential",
      TRUE ~ Decay_Type
    )
  )


# Plot
viral_decay_plot <- ggplot(viral_decay_long %>% dplyr::filter(Population == 'c_Viruses'), aes(x = combi_tag, y = Value, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = 'black') +
  facet_wrap(Decay_Type ~ Decay_Metric, scales = "free") +
  labs(x = "Stations", y = "Value", title = "Decay Rate and Percent Decay") +
  theme_test() +
  #scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
viral_decay_plot
ggsave(plot = viral_decay_plot, filename = "./figures/viral_decay_plot.svg", width = 12, height = 8, dpi = 400)

