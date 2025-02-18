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
  "Chlorophyll" = "Chlorophyll~(mu*g~L^{-1})",
  "HNA" = "HNA~(10^6~cells~mL^{-1})",
  "LNA" = "LNA~(10^6~cells~mL^{-1})",
  "Cyanobacteria" = "Cyanobacteria~(10^3~cells~mL^{-1})",
  "V1" = "V1~(10^6~VLPs~mL^{-1})",
  "V2" = "V2~(10^6~VLPs~mL^{-1})",
  "V3" = "V3~(10^6~VLPs~mL^{-1})",
  "Corrected_VP_Lytic" = "Lytic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "Corrected_VP_Lysogenic" = "Lysogenic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "Corrected_VP_Lytic_per_bacteria" = "Normalized~lytic~production~rate~(h^{-1})",
  "Corrected_VP_Lysogenic_per_bacteria" = "Normalized~lysogenic~production~rate~(~h^{-1})",
  "decay_rate_linear" = "Viral~decay~rate~(10^3~cells~mL^{-1}~h^{-1})",
  "percent_bacterial_loss_day_burst_50" = "Bacteria~lysed~('%'~cells~d^{-1})",
  "percent_lysogeny_day_burst_50" = "Lysogeny~('%'~cells~d^{-1})",
  "percent_decay_day_linear" = "Viral~decay~('%'~VLPs~d^{-1})"
)

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


pc_boxplots_3d <- ggplot(physicochemical_df_3d,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `physicochemical parameters`, scales = "free_y", ncol = 3, , 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = custom_color_palette_cruise) +
  labs(
    fill = "Seasons"
  ) +
  theme_test(base_size = 15) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

pc_boxplots_3d

ggsave(plot = pc_boxplots_3d, filename = "./figures/boxplots_physicochemical_3depths.svg", dpi = 800, width = 12, height = 7)

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


bio_boxplots_3d <- ggplot(biological_df_3d,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `biological parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = custom_color_palette_cruise) +
  labs(
    fill = "Seasons"
  ) +
  theme_test(base_size = 15) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

bio_boxplots_3d

ggsave(plot = bio_boxplots_3d, filename = "./figures/boxplots_biological_3depths.svg", dpi = 800, width = 12, height = 7)



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


pc_boxplots_7m <- ggplot(physicochemical_df_7m,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `physicochemical parameters`, scales = "free_y",  ncol = 3, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = custom_color_palette_cruise) +
  labs(
    fill = "Seasons"
  ) +
  theme_test(base_size = 15) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

pc_boxplots_7m

ggsave(plot = pc_boxplots_7m, filename = "./figures/boxplots_physicochemical_7m.svg", dpi = 800, width = 12, height = 7)

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


bio_boxplots_7m <- ggplot(biological_df_7m,
                          aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `biological parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = custom_color_palette_cruise) +
  labs(
    fill = "Seasons"
  ) +
  theme_test(base_size = 15) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

bio_boxplots_7m

ggsave(plot = bio_boxplots_7m, filename = "./figures/boxplots_biological_7m.svg", dpi = 800, width = 12, height = 7)


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


vp_boxplots_7m <- ggplot(vp_df_7m,
                          aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `vp parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(variable_labels, label_parsed) #, strip.position = "bottom"
  ) + 
  scale_fill_manual(values = custom_color_palette_cruise) +
  labs(
    fill = "Seasons"
  ) +
  theme_test(base_size = 15) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") 

vp_boxplots_7m

ggsave(plot = vp_boxplots_7m, filename = "./figures/boxplots_vp_7m.svg", dpi = 800, width = 12, height = 7)

