library(tidyverse)
library(ggsci)


setwd("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production")
#Import FCM count csv
counts_per_mL<- read.csv("./PE_Cruises/results/PE_Cruises_per_mL.csv")
str(counts_per_mL)
counts_0.22<- counts_per_mL %>%
  dplyr::filter(Sample_Type == '0.22',
                Staining_Protocol == 'Viruses') %>%
  select(Location, Station_Number, Depth, Sample_Type, Timepoint, Replicate, c_Viruses, c_V1, c_V2, c_V3) %>%
  pivot_longer(cols = c(c_Viruses, c_V1, c_V2, c_V3),
               names_to = "Population",
               values_to = "Counts") %>%
  group_by(Location, Station_Number, Depth, Sample_Type, Timepoint, Population) %>%
  summarise(Mean_Counts = mean(Counts, na.rm = TRUE), .groups = 'drop')

# Calculating decay rates

viral_decay <- counts_0.22 %>%
  pivot_wider(names_from = "Timepoint",
              values_from = "Mean_Counts") %>%
  mutate(decay_rate = (`24` - `0`)/24)

# Adding orignal viral abundance to calulate percent viruses lost per day
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
  mutate(percent_decay_day = (decay_rate/Original_Abundance)*100*24) %>%
  select(-`0`, -`24`, -Original_Abundance)


write.csv(viral_decay_percent, "./results/viral_decay.csv", row.names = F, quote = F)


viral_decay <- viral_decay %>%
  mutate(combi_tag = paste(Location, Station_Number, sep = "_"))

# Pivot longer for plotting
viral_decay_long <- viral_decay %>%
  pivot_longer(cols = c(decay_rate, percent_decay_day), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(Population = factor(Population, levels = c("c_Viruses", "c_V1", "c_V2", "c_V3")))

# Plot
ggplot(viral_decay_long, aes(x = combi_tag, y = Value, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Metric, scales = "free_y", ncol =1) +
  labs(x = "Combi Tag", y = "Value", title = "Decay Rate and Percent Decay by Population") +
  theme_minimal() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
