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

pe486_df <- pe_df %>%
  dplyr::filter(Location == "PE486")

# Assigning Dogger bank group
pe486_3d_df <- pe486_df %>%
  mutate(Dogger = case_when(
    Station_Number %in% c(9.0, 12.1, 12.2) ~ "Dogger's Bank",
    TRUE ~ "Non-Dogger's Bank"
  ))


# 2.0 Assigning vectors for parameters ####

physicochemical_params <- c("Temperature", "Salinity", "Density", "Conductivity", 
                            "Turbidity", "Nitrate", "Phosphate", "Silicate", "Max_Depth")

biological_params <- c("Oxygen", "Fluorescence",
                       "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
                       "Total_Viruses", "V1", "V2", "V3",
                       "VBR", "HNALNA")

vp_params <- c("decay_rate_linear",
               "percent_decay_day_linear",
               "Corrected_VP_Lytic",
               "Corrected_VP_Lysogenic",
               # "Corrected_VP_SE_Lytic",
               # "Corrected_VP_SE_Lysogenic",
               "percent_bacterial_loss_day_burst_50",
               "percent_lysogeny_day_burst_50")

sample_params <- c("Location", "Station_Number", "Depth", 
                   "sample_tag", "Latitude", "Longitude",
                   "Season")


# Defining labels using parse-able expressions
variable_labels <- c(
  "Temperature" = "Temperature~(degree*C)",
  "Salinity" = "Salinity~(PSU)",
  "Turbidity" = "Turbidity~(NTU)",
  "Conductivity" = "Conductivity~(mS/cm)",
  "Density" = "Density~(kg/m^3)", 
  "Max_Depth" = "Max~Depth~(m)",
  "Nitrate" = "Nitrate~(mu*M)",
  "Phosphate" = "Phosphate~(mu*M)",
  "Silicate" = "Silicate~(mu*M)",
  "Total_Bacteria" = "Bacteria~(10^6~cells~mL^{-1})",
  "Total_Viruses" = "Viruses~(10^6~VLPs~mL^{-1})",
  "VBR" = "VBR",
  "HNALNA" = "HNALNA",
  "Oxygen" = "Oxygen~(mu*M)",  
  "Fluorescence" = "Fluorescence~(mu*g~L^{-1})",
  "HNA" = "HNA~(10^6~cells~mL^{-1})",
  "LNA" = "LNA~(10^6~cells~mL^{-1})",
  "Cyanobacteria" = "Cyanobacteria~(10^3~cells~mL^{-1})",
  "V1" = "V1~(10^6~VLPs~mL^{-1})",
  "V2" = "V2~(10^6~VLPs~mL^{-1})",
  "V3" = "V3~(10^6~VLPs~mL^{-1})",
  "Corrected_VP_Lytic" = "Lytic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "Corrected_VP_Lysogenic" = "Lysogenic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "decay_rate_linear" = "Viral~decay~rate~(10^3~cells~mL^{-1}~h^{-1})",
  "percent_bacterial_loss_day_burst_50" = "Bacteria~lysed~('%'~cells~d^{-1})",
  "percent_lysogeny_day_burst_50" = "Lysogeny~('%'~cells~d^{-1})",
  "percent_decay_day_linear" = "Viral~decay~('%'~VLPs~d^{-1})"
)
# 3.0 Boxplots for 3 depths####


doggers_3depths_df <- pe486_3d_df %>%
  dplyr::select(c(all_of(c(physicochemical_params, sample_params, biological_params)), "Dogger")) %>%
  pivot_longer(cols = all_of(c(physicochemical_params, biological_params)),
               names_to = "parameters",
               values_to = "Value") %>%
  dplyr::mutate(`parameters` = factor(`parameters`, levels = c("Temperature", "Salinity", "Conductivity", "Density", 
                                                               "Turbidity", "Nitrate", "Phosphate", "Silicate",
                                                               "Total_Bacteria", "HNA", "LNA", "Cyanobacteria", 
                                                               "Total_Viruses","V1", "V2", "V3",
                                                               "VBR", "HNALNA", "Oxygen", "Fluorescence", "Max_Depth"
                                                               )))

doggers_3d_boxplots <- ggplot(doggers_3depths_df,
                         aes(x = Dogger, y = Value, fill = Dogger)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `parameters`, scales = "free_y", nrow = 4, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = c(`Dogger's Bank` = "#2a9d8f", `Non-Dogger's Bank` = "#37474f")) +
  labs(
    fill = "Dogger"
  ) +
  theme_test(base_size = 20) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

doggers_3d_boxplots

ggsave(plot = doggers_3d_boxplots, filename = "./figures/boxplots_doggers_bank_3depths.svg", dpi = 800, width = 30, height = 12)

# Perform Kruskal-Wallis test for all parameters

kw_results_doggers_3d <- doggers_3depths_df %>%
  group_by(`parameters`) %>%
  kruskal_test(Value ~ Dogger) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_doggers_3d, n = 30)



# 4.0 Boxplots for 7m####


doggers_7m_df <- pe486_3d_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(c(all_of(c(physicochemical_params, sample_params, biological_params, vp_params)), "Dogger")) %>%
  pivot_longer(cols = all_of(c(physicochemical_params, biological_params, vp_params)),
               names_to = "parameters",
               values_to = "Value") %>%
  dplyr::mutate(`parameters` = factor(`parameters`, levels = c("Temperature", "Salinity", "Conductivity", "Density", 
                                                               "Turbidity", "Nitrate", "Phosphate", "Silicate",
                                                               "Total_Bacteria", "HNA", "LNA", "Cyanobacteria", 
                                                               "Total_Viruses","V1", "V2", "V3",
                                                               "VBR", "HNALNA", "Oxygen", "Fluorescence",
                                                              "Corrected_VP_Lytic", "Corrected_VP_Lysogenic", "decay_rate_linear", "Max_Depth", 
                                                               "percent_bacterial_loss_day_burst_50","percent_lysogeny_day_burst_50",
                                                               "percent_decay_day_linear")))


doggers_7m_boxplots <- ggplot(doggers_7m_df,
                              aes(x = Dogger, y = Value, fill = Dogger)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `parameters`, scales = "free_y", ncol = 4, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = c(`Dogger's Bank` = "#2a9d8f", `Non-Dogger's Bank` = "#37474f")) +
  labs(
    fill = "Dogger"
  ) +
  theme_test(base_size = 20) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

doggers_7m_boxplots

ggsave(plot = doggers_7m_boxplots, filename = "./figures/boxplots_doggers_bank_7m.svg", dpi = 800, width = 20, height = 15)

# Perform Kruskal-Wallis test for all parameters

kw_results_doggers_7m <- doggers_7m_df %>%
  group_by(`parameters`) %>%
  kruskal_test(Value ~ Dogger) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_doggers_7m, n = 30)
