#AIM: To check how viralprod VIPCAL and VIPCAL-SE compare with linear regression approach

# 0.0 Setting up ####
source("scripts/0_source.R")


# 1.0 Importing viral production rate data frames ####

viralprod <- read.csv("results/viral_production_analyses/PE_Cruises_viral_production/vp_results_ALL.csv")

#Linear regression output file is saved in metadata

lm <- readxl::read_excel("metadata/PE477_PE486_VP_LM_HMS.xlsx", sheet = 'final_output')

# 2.0 Wrangling and combining the two data frames ####


viralprod_transformed <- viralprod %>% 
  dplyr::select(-c("abs_VP", "VP_SE", "VP_R_Squared")) %>%
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
lm_vp <- lm %>%
  select(Location, Station_Number, Depth, corrected_Lytic, corrected_Lysogenic) %>%
  pivot_longer(
    cols = c(corrected_Lytic, corrected_Lysogenic),
    names_to = "Treatment",
    values_to = "VP"
  ) %>%
  mutate(
    Treatment = recode(Treatment,
                       corrected_Lytic = "Lytic",
                       corrected_Lysogenic = "Lysogenic")
  )

# Splitting  data for Time_Ranges (Lytic_Time_Range and Lysogenic_Time_Range)
lm_time <- lm %>%
  select(Location, Station_Number, Depth, Lytic_Time_Range, Lysogenic_Time_Range) %>%
  pivot_longer(
    cols = c(Lytic_Time_Range, Lysogenic_Time_Range),
    names_to = "Treatment",
    values_to = "Time_Range"
  ) %>%
  mutate(
    Treatment = recode(Treatment,
                       Lytic_Time_Range = "Lytic",
                       Lysogenic_Time_Range = "Lysogenic"),
    VP = if_else(column_name < 0, 0, VP)
  )

# Joining dataframes 
lm_transformed <- lm_vp %>%
  dplyr::left_join(lm_time, by = c("Location", "Station_Number", "Depth", "Treatment")) %>%
  dplyr::mutate(VP_Method = "LR",
                Population = "c_Viruses")

# Combinign viralprod and LM

vp <- rbind(viralprod_transformed, lm_transformed)

# Adding viral production assay bacterial efficiency information 
vp <- vp %>%
  left_join(
    lm %>% select(Location, Station_Number, Depth, bac_efficiency), 
    by = c("Location", "Station_Number", "Depth")
  ) %>%
  mutate(
    VP_Method = factor(VP_Method, levels = c("LR", "VIPCAL", "VIPCAL-SE")),
    Treatment = factor(Treatment, levels = c("Lytic", "Lysogenic"))
  )

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

# 4.0 Hypothesis testing ####

vp_lytic <- vp %>% dplyr::fiter(Treatment == 'Lytic')
vp_lysogenic <- vp %>% dplyr::fiter(Treatment == 'Lysogenic')
# Testing for normality

shapiro.test(vp$VP) #not normal
shapiro.test(vp_lytic$VP) #not normal
shapiro.test(vp_lysogenic$VP) #not normal

# Lytic production
kruskal.test(VP ~ VP_Method, data = vp_lytic) 
#Kruskal-Wallis chi-squared = 6.9802, df = 2, p-value = 0.0305
# significant difference
pairwise.wilcox.test(vp_lytic$VP, vp_lytic$VP_Method, p.adjust.method = "bonferroni")
# only between LR and VIPCAL. VIPCAL-SE is not different from either
kruskal.test(VP ~ VP_Method, data = vp_lysogenic)
# Kruskal-Wallis chi-squared = 3.1328, df = 2, p-value = 0.2088
# no significant difference
pairwise.wilcox.test(vp_lysogenic$VP, vp_lysogenic$VP_Method, p.adjust.method = "bonferroni")
