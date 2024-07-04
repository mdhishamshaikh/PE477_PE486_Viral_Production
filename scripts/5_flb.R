library(tidyverse)
library(ggsci)
library(readxl)
library(lubridate)


# Importing the dataframe
flb_df <- read_excel("./Data/flb_pe477/FLB_PE477.xlsx", 
                     sheet = "Sheet1")


# Selcting for counts for PE477
flb_df <- flb_df %>%
  filter(Sample_type != 'Stock',
         Location == 'PE477')


flb_df_T0 <- flb_df %>% filter(Time == 0) %>% select(Day, Location, Sample_type, Replicate,cells_per_mL)
flb_df_T24 <- flb_df %>% filter(Time == 24) %>% select(Day, Location, Sample_type, Replicate, cells_per_mL)

flb_df_merged <- full_join(flb_df_T0, flb_df_T24, by = c("Day", "Location", "Replicate", "Sample_type"), suffix = c("_T0", "_T24"))

flb_df_merged <- flb_df_merged %>%
  mutate(cells_per_mL_diff = cells_per_mL_T24 - cells_per_mL_T0,
         sample_rep = paste0(Sample_type, "_", Replicate))


ggplot(flb_df_merged, aes(x = Sample_type, y = cells_per_mL_diff, fill = sample_rep)) +
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
  filter(Sample_type != "Control") %>%
  group_by(Day, Sample_type) %>%
  summarise(
    mean_diff = mean(cells_per_mL_diff, na.rm = TRUE),
    sd_diff = sd(cells_per_mL_diff, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate mean and SD for control samples
control_stats <- flb_df_merged %>%
  filter(Sample_type == "Control") %>%
  group_by(Day) %>%
  summarise(
    mean_diff = mean(cells_per_mL_diff, na.rm = TRUE),
    sd_diff = sd(cells_per_mL_diff, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(Sample_type = "Control")

# Combine FLB and control stats
stats_combined <- bind_rows(flb_stats, control_stats)


ggplot(stats_combined, aes(x = interaction(Sample_type, Day), y = mean_diff, fill = Sample_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_npg() +
  labs(title = "Average Cells per mL Difference (T24 - T0) by Sample Type and Day",
       x = "Sample Type and Day",
       y = "Average Cells per mL Difference",
       fill = "Sample Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

control_stats <- control_stats%>%
  rename(control_mean_diff = mean_diff)

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
