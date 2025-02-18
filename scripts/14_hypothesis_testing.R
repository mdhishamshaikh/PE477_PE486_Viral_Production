# AIM: To create boxplots #####

# Setting up ####

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv") %>%
  dplyr::mutate(Total_Bacteria = Total_Bacteria/1e+6,
                Total_Viruses = Total_Viruses/1e+6,
                HNA = HNA/1e+6,
                LNA = LNA/1e+6,
                V1 = V1/1e+6,
                V2 = V2/1e+6,
                V3 = V3/1e+6,
                Cyanobacteria = Cyanobacteria/1e+3,
                decay_rate_linear  = decay_rate_linear /1e+3,
                Corrected_VP_Lytic  = Corrected_VP_Lytic/1e+5,
                Corrected_VP_Lysogenic = Corrected_VP_Lysogenic/1e+5
  )


# 2.0 Assigning vectors for parameters ####

physicochemical_params <- c("Temperature", "Salinity", #"Density", "Conductivity", 
                            "Turbidity", "Nitrate", "Phosphate", "Silicate")

biological_params <- c(#"Oxygen", 
  "Chlorophyll",
  "Total_Bacteria", #"HNA", "LNA", 
  "Cyanobacteria",
  "Total_Viruses", #"V1", "V2", "V3",
  "VBR")

vp_params <- c("decay_rate_linear",
               "percent_decay_day_linear",
               "Corrected_VP_Lytic",
               "Corrected_VP_Lysogenic",
               # "Corrected_VP_SE_Lytic",
               # "Corrected_VP_SE_Lysogenic",
               "Corrected_VP_Lytic_per_bacteria",
               "Corrected_VP_Lysogenic_per_bacteria",
               "percent_bacterial_loss_day_burst_50",
               "percent_lysogeny_day_burst_50")

sample_params <- c("Location", "Station_Number", "Depth", 
                   "sample_tag", "Latitude", "Longitude",
                   "Season")
# 3.0 3 depths ####

# 3.1 Physicochemical parameters ####

physicochemical_df_3d <- pe_df %>%
  dplyr::select(all_of(c(physicochemical_params, sample_params))) %>%
  pivot_longer(cols = all_of(physicochemical_params),
               names_to = "physicochemical parameters",
               values_to = "Value") %>%
  dplyr::mutate(`physicochemical parameters` = factor(`physicochemical parameters`, levels = c("Temperature", "Salinity", 
                                                                                               "Turbidity", "Nitrate", "Phosphate", "Silicate"))) %>%
  na.omit()


# Performing Wilcoxon test for each physicochemical parameter
wilcox_results_pc_3d <- physicochemical_df_3d %>%
  group_by(`physicochemical parameters`) %>%
  wilcox_test(Value ~ Season, exact = TRUE) %>%  #`exact = TRUE` for small sample sizes
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns" 
  )) %>%
  select(`physicochemical parameters`, statistic, p, Significance) %>%
  rename(
    Parameter = `physicochemical parameters`,
    W_Statistic = statistic,
    P_Value = p
  )

print(wilcox_results_pc_3d)

# 3.2 Biological parameters ####

biological_df_3d <- pe_df %>%
  dplyr::select(all_of(c(biological_params, sample_params))) %>%
  pivot_longer(cols = all_of(biological_params),
               names_to = "biological parameters",
               values_to = "Value") %>%
  dplyr::mutate(`biological parameters` = factor(`biological parameters`, levels = c("Total_Bacteria", "Total_Viruses", "VBR", 
                                                                                     "Chlorophyll",
                                                                                     "Cyanobacteria"))) %>%
  na.omit()

# Performing Wilcoxon test for each biological parameter
wilcox_results_bio_3d <- biological_df_3d %>%
  group_by(`biological parameters`) %>%
  wilcox_test(Value ~ Season, exact = TRUE) %>%  #`exact = TRUE` for small sample sizes
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns" 
  )) %>%
  select(`biological parameters`, statistic, p, Significance) %>%
  rename(
    Parameter = `biological parameters`,
    W_Statistic = statistic,
    P_Value = p
  )

print(wilcox_results_bio_3d)

# 4.0 7 m ####

# 4.1 Physicochemical parameters ####

physicochemical_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(physicochemical_params, sample_params))) %>%
  pivot_longer(cols = all_of(physicochemical_params),
               names_to = "physicochemical parameters",
               values_to = "Value") %>%
  dplyr::mutate(`physicochemical parameters` = factor(`physicochemical parameters`, levels = c("Temperature", "Salinity", 
                                                                                               "Turbidity", "Nitrate", "Phosphate", "Silicate"))) %>%
  na.omit()


# Performing Wilcoxon test for each physicochemical parameter
wilcox_results_pc_7m <- physicochemical_df_7m %>%
  group_by(`physicochemical parameters`) %>%
  wilcox_test(Value ~ Season, exact = TRUE) %>%  #`exact = TRUE` for small sample sizes
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns" 
  )) %>%
  select(`physicochemical parameters`, statistic, p, Significance) %>%
  rename(
    Parameter = `physicochemical parameters`,
    W_Statistic = statistic,
    P_Value = p
  )

print(wilcox_results_pc_7m)

# 4.2 Biological parameters ####

biological_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(biological_params, sample_params))) %>%
  pivot_longer(cols = all_of(biological_params),
               names_to = "biological parameters",
               values_to = "Value") %>%
  dplyr::mutate(`biological parameters` = factor(`biological parameters`, levels = c("Total_Bacteria", "Total_Viruses", "VBR", 
                                                                                     "Chlorophyll",
                                                                                     "Cyanobacteria"))) %>%
  na.omit()


# Performing Wilcoxon test for each biological parameter
wilcox_results_bio_7m <- biological_df_7m %>%
  group_by(`biological parameters`) %>%
  wilcox_test(Value ~ Season, exact = TRUE) %>%  #`exact = TRUE` for small sample sizes
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns" 
  )) %>%
  select(`biological parameters`, statistic, p, Significance) %>%
  rename(
    Parameter = `biological parameters`,
    W_Statistic = statistic,
    P_Value = p
  )

print(wilcox_results_bio_7m)

# 4.3 vp parameters ####

vp_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(vp_params, sample_params))) %>%
  pivot_longer(cols = all_of(vp_params),
               names_to = "vp parameters",
               values_to = "Value") %>%
  dplyr::mutate(`vp parameters` = factor(`vp parameters`, levels = c( "Corrected_VP_Lytic", "percent_bacterial_loss_day_burst_50", "Corrected_VP_Lytic_per_bacteria", 
                                                                      "Corrected_VP_Lysogenic", "percent_lysogeny_day_burst_50", "Corrected_VP_Lysogenic_per_bacteria", 
                                                                      "decay_rate_linear",  "percent_decay_day_linear" )),
                Value = ifelse(Value == 0, NA, Value)) %>%
  na.omit()


# Performing Wilcoxon test for each VP parameter
wilcox_results_vp_7m <- vp_df_7m %>%
  group_by(`vp parameters`) %>%
  wilcox_test(Value ~ Season, exact = TRUE) %>%  #`exact = TRUE` for small sample sizes
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns" 
  )) %>%
  select(`vp parameters`, statistic, p, Significance) %>%
  rename(
    Parameter = `vp parameters`,
    W_Statistic = statistic,
    P_Value = p
  )

print(wilcox_results_vp_7m)


# 5.0 Combining Wilcoxon results ####

wilcoxon_results <- ls(pattern = "^wilcox")
combined_wilcoxon_results <- bind_rows(
  lapply(wilcoxon_results, function(var) {
    wilcoxon_df <- get(var)  
    wilcoxon_df <- wilcoxon_df %>% dplyr::mutate(`Tested For` = var)  
    return(wilcoxon_df)
  })
)

print(combined_wilcoxon_results, n = 50)
