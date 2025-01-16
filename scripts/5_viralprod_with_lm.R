#AIM: To check how viralprod VIPCAL and VIPCAL-SE compare with linear regression approach

# 0.0 Setting up ####
source("scripts/0_source.R")


# 1.0 Importing viral production rate data frames ####

viralprod <- read.csv("results/viral_production_analyses/PE_Cruises_viral_production/vp_results_ALL.csv")

#Linear regression output file is saved in metadata

lr <- readxl::read_excel("metadata/PE477_PE486_VP_LM_HMS.xlsx", sheet = 'final_output')

# 2.0 Wrangling and combining the two data frames ####


viralprod_transformed <- viralprod %>% 
  dplyr::select(-c("VP_R_Squared")) %>%
  dplyr::filter(Population == "c_Viruses",
                Sample_Type %in% c('VP', 'Diff'),
                Time_Range == "T0_T24", # The NGTE is 24 hours for all stations (low bacterial growth)
                VP_Method %in% c("VPCL_AR_DIFF", "VPCL_AR_DIFF_SE"))  %>%
  dplyr::mutate(
    Depth = ifelse(Depth == 1, 7, Depth), 
    Sample_Type = dplyr::recode(Sample_Type, 
                                VP = "Lytic", 
                                Diff = "Lysogenic"), 
    VP_Method = dplyr::recode(VP_Method, 
                              VPCL_AR_DIFF = "VIPCAL", 
                              VPCL_AR_DIFF_SE = "VIPCAL-SE")
  ) %>%
  dplyr::rename(Treatment = Sample_Type)

# Splitting  data for VP values (corrected_Lytic and corrected_Lysogenic)
lr_vp <- lr %>%
  select(Location, Station_Number, Depth, corrected_Lytic, corrected_Lysogenic) %>%
  pivot_longer(
    cols = c(corrected_Lytic, corrected_Lysogenic),
    names_to = "Treatment",
    values_to = "VP"
  ) %>%
  mutate(
    Treatment = dplyr::recode(Treatment,
                       "corrected_Lytic" = "Lytic",
                       "corrected_Lysogenic" = "Lysogenic")
  )

# Splitting  data for Time_Ranges (Lytic_Time_Range and Lysogenic_Time_Range)
lr_time <- lr %>%
  select(Location, Station_Number, Depth, Lytic_Time_Range, Lysogenic_Time_Range) %>%
  pivot_longer(
    cols = c(Lytic_Time_Range, Lysogenic_Time_Range),
    names_to = "Treatment",
    values_to = "Time_Range"
  ) %>%
  mutate(
    Treatment = dplyr::recode(Treatment,
                       Lytic_Time_Range = "Lytic",
                       Lysogenic_Time_Range = "Lysogenic"),
    VP_SE = NA,
    abs_VP = NA)

# Joining dataframes 
lr_transformed <- lr_vp %>%
  dplyr::left_join(lr_time, by = c("Location", "Station_Number", "Depth", "Treatment")) %>%
  dplyr::mutate(VP_Method = "LR",
                Population = "c_Viruses",
                VP = if_else(VP < 0, 0, VP))

# Combinign viralprod and LM

vp <- rbind(viralprod_transformed, lr_transformed)

# Adding viral production assay bacterial efficiency information 
vp <- vp %>%
  left_join(
    lm %>% select(Location, Station_Number, Depth, bac_efficiency), 
    by = c("Location", "Station_Number", "Depth")
  ) %>%
  mutate(
    VP_Method = factor(VP_Method, levels = c("LR", "VIPCAL", "VIPCAL-SE")),
    Treatment = factor(Treatment, levels = c("Lytic", "Lysogenic")),
    Location_Station = paste(Location, Station_Number, sep = "_")
  )

write.csv(vp, "./results/viral_production_analyses/viralproduction_viralprod_lr.csv", row.names = F)

# 3.0 visualization ####
#to see if there is a difference between lytic and lysogenic viral production between the methods

vp_methods_plot <- ggplot(vp, aes(x = VP_Method, y = VP, fill = VP_Method)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Location), width = 0.2, alpha = 0.7) + 
  facet_grid(. ~Treatment )+
  labs(
    title = "Viral production across methods",
    x = "VP Method",
    y = "Viral Production (VP)"
  ) +
  theme_test() +
  scale_fill_npg() +
  scale_color_manual(values = c("PE477" = "#0c1844", "PE486" = "#850000"))
vp_methods_plot

ggsave(plot = vp_methods_plot, filename = "./figures/vp_methods.svg", width = 6, height = 4, dpi = 400)


# Grouped bar plot
vp_grouped_barplot <- ggplot(vp, aes(x = Location_Station, y = VP, fill = VP_Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7#, color = "black"
           ) +
  facet_wrap(~ Treatment, scales = "free_y") +  # Facet by Treatment
  labs(
    title = "Viral Production by Station and Method",
    x = "Stations",
    y =  expression("Viral production rate" ~ (x ~ 10^6 ~ "VLPs" ~ mL^-1 ~ h^-1)),
    fill = "VP Method"
  ) +
  theme_test(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  scale_fill_manual(values = c("LR" = "#DC2828", "VIPCAL" = "#003049", "VIPCAL-SE" ="#2A9D8F", "VIPCAL-SE-GTE" ="#F6CD61"))

vp_grouped_barplot
ggsave(plot = vp_grouped_barplot, filename = "./figures/vp_methods_by_stations.svg", width = 12, height = 4, dpi = 400)


# 4.0 Hypothesis testing ####

vp_lytic <- vp %>% dplyr::filter(Treatment == 'Lytic')
vp_lysogenic <- vp %>% dplyr::filter(Treatment == 'Lysogenic')
# Testing for normality

shapiro.test(vp$VP) #not normal
shapiro.test(vp_lytic$VP) #not normal
shapiro.test(vp_lysogenic$VP) #not normal

# Lytic production
kruskal.test(VP ~ VP_Method, data = vp_lytic) 
#Kruskal-Wallis chi-squared = 7.1217, df = 2, p-value = 0.02841
# significant difference
pairwise.wilcox.test(vp_lytic$VP, vp_lytic$VP_Method, p.adjust.method = "bonferroni")
# only between LR and VIPCAL. VIPCAL-SE is not different from either
kruskal.test(VP ~ VP_Method, data = vp_lysogenic)
# Kruskal-Wallis chi-squared = 2.4904, df = 2, p-value = 0.2879
# no significant difference
pairwise.wilcox.test(vp_lysogenic$VP, vp_lysogenic$VP_Method, p.adjust.method = "bonferroni")
