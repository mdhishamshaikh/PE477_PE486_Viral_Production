library(tidyverse)
library(ggsci)
library(readxl)
library(lubridate)


# Importing the dataframe
flb_df <- read_excel("./Data/flb_pe477/FLB_PE477.xlsx", 
                     sheet = "Sheet1")


# Selcting for counts for PE477
flb_df <- flb_df %>%
  dplyr::rename(Sample_Type = Sample_type) %>%
  filter(Sample_Type != 'Stock',
         Location == 'PE477')


flb_df_T0 <- flb_df %>% filter(Time == 0) %>% select(Day, Location, Sample_Type, Replicate,cells_per_mL)
flb_df_T24 <- flb_df %>% filter(Time == 24) %>% select(Day, Location, Sample_Type, Replicate, cells_per_mL)

flb_df_merged <- full_join(flb_df_T0, flb_df_T24, by = c("Day", "Location", "Replicate", "Sample_Type"), suffix = c("_T0", "_T24"))

flb_df_merged <- flb_df_merged %>%
  mutate(
    cells_per_mL_diff = cells_per_mL_T24 - cells_per_mL_T0,
    sample_rep = paste0(Sample_Type, "_", Replicate),
    exp_loss_rate = (log(cells_per_mL_T24) - log(cells_per_mL_T0)) / 24  # Assuming time step is 24 hours
  )


ggplot(flb_df_merged, aes(x = Sample_Type, y = cells_per_mL_diff, fill = sample_rep)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Day) +
  scale_fill_npg() +
  labs(title = "Difference in Cells per mL (T24 - T0) by Sample Type",
       x = "Sample Type",
       y = "Cells per mL Difference",
       fill = "Sample Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Calculate mean and SD for FLB sample types
flb_stats <- flb_df_merged %>%
  dplyr::filter(Sample_Type != "Control") %>%
  group_by(Day) %>%
  dplyr::summarise(
    mean_diff = mean(cells_per_mL_diff, na.rm = TRUE),
    se_diff = plotrix::std.error(cells_per_mL_diff, na.rm = TRUE),
    exp_loss_rate_mean = mean(exp_loss_rate, na.rm = T) 
  ) %>%
  ungroup() %>%
  mutate(Sample_Type = "FLB")

# Calculate mean and SD for control samples
control_stats <- flb_df_merged %>%
  filter(Sample_Type == "Control") %>%
  group_by(Day) %>%
  dplyr::summarise(
    mean_diff = mean(cells_per_mL_diff, na.rm = TRUE),
    se_diff = plotrix::std.error(cells_per_mL_diff, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(Sample_Type = "Control")

# Combine FLB and control stats
stats_combined <- bind_rows(flb_stats, control_stats) %>%
  dplyr::rename("Station_Number" = "Day") %>%
  mutate(Location = "PE477",
         Depth = 1)


# Adding orignal bacterial abundance to calulate percent bacteria lost per day
abundance_ts <- read.csv("./Data/nutrients_ts_abundance_coordinates/nutrients_ts_abundance_coordinates.csv")
abundance <- abundance_ts %>%
  dplyr::select(Location, Station_Number, Depth, Total_Bacteria) %>%
  dplyr::rename(c_Bacteria = Total_Bacteria)

flb_grazing_abundance<- stats_combined %>%
  left_join(abundance, by = c("Location", "Station_Number", "Depth")) %>%
  dplyr::rename(flb_loss_mean_rate = mean_diff,
         flb_loss_se = se_diff)

flb_grazing_percent <- flb_grazing_abundance %>%
  mutate(percent_flb_grazed_day = (flb_loss_mean_rate/c_Bacteria)*100*24) %>%
  dplyr::select(-c_Bacteria)


write.csv(flb_grazing_percent, "./results/flb_grazing.csv", row.names = F, quote = F)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   




# The plots need adjustment wrt to their names

flb_control_stats_plot<- ggplot(stats_combined, aes(x = interaction(Sample_Type, Station_Number), y = mean_diff, fill = Sample_Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_npg() +
  labs(title = "FLB grazed by heteronanoflagellates",
       x = "Sample Type and Day",
       y = "Average Cells per mL Difference",
       fill = "Sample Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
flb_control_stats_plot
ggsave(flb_control_stats_plot, path = "./results/figures/", filename = "flb_control_PE477.svg", dpi = 800, height = 3, width = 4.5)

control_stats <- control_stats%>%
  dplyr::rename(control_mean_diff = mean_diff)

stats_combined <- flb_stats %>%
  merge(control_stats, by = "Day") %>%
  mutate(diff_from_control = mean_diff - control_mean_diff)

ggplot(stats_combined, aes(x = Day, y = diff_from_control)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~ Day) +
  scale_fill_npg() +
  labs(title = "Difference in Cells per mL Difference (FLB - Control) by Day",
       x = "Sample Type",
       y = "Difference from Control") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
