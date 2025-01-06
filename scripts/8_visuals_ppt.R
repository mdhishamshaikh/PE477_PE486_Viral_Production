library(tidyverse)
library(ggsci)


vp_pe_df <- read.csv("./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_analyzed_pe_cruises_combined.txt", sep = "\t")


#### 1.0 Bacterial abundance ####

bac_long <- vp_pe_df %>%
  select(Cruise, Station, Total_Bacteria, HNA, LNA) %>%
  pivot_longer(cols = c(Total_Bacteria, HNA, LNA),
               names_to = "Population",
               values_to = "Count") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),  
         Population = factor(Population, levels = c("Total_Bacteria", "HNA", "LNA"))) %>%
  mutate(Count = Count / 1e6)


bac_count_plot<- ggplot(bac_long, aes(x = Station, y = Count, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  labs(x = "Station",
       y = "Bacterial abundance (x 10^6 cells mL-1)",
       fill = "Population") +
  scale_fill_manual(values = c("Total_Bacteria" = "#850000", "HNA" = "#B80000", "LNA" = "#FF6969")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside") +
  guides(x = guide_axis(position = "bottom")) +  
  scale_x_continuous(breaks = bac_long$Station, labels = bac_long$Station)

bac_count_plot

ggsave(bac_count_plot, path = "./results/figures/", filename = "bacterial_abundance_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)

#### 2.0 Viral abundance ####

vir_long <- vp_pe_df %>%
  select(Cruise, Station, Total_Viruses, V1, V2, V3) %>%
  pivot_longer(cols = c(Total_Viruses, V1, V2, V3),
               names_to = "Population",
               values_to = "Count") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),  
         Population = factor(Population, levels = c("Total_Viruses", "V1", "V2", "V3"))) %>%
  mutate(Count = Count / 1e6)


vir_count_plot<- ggplot(vir_long, aes(x = Station, y = Count, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  labs(x = "Station",
       y = "Viral abundance (x 10^6 cells mL-1)",
       fill = "Population") +
  scale_fill_manual(values = c("Total_Viruses" = "#0c1844", "V3" = "#4d778b", "V2" = "#6baed6", "V1" = "#474F7A")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside") +
  guides(x = guide_axis(position = "bottom")) +  
  scale_x_continuous(breaks = vir_long$Station, labels = vir_long$Station)

vir_count_plot

ggsave(vir_count_plot, path = "./results/figures/", filename = "viral_abundance_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)


#### 3.0 Abundance ratio ####


vbr_df <- vp_pe_df %>%
  select(Cruise, Station, VBR) %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486"))) 

vbr_plot <- ggplot(vbr_df, aes(x = Station, y = VBR)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#17153B") +
  labs(x = "Station",
       y = "VBR",
       title = "VBR") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside",
        legend.position = "none") +
  guides(x = guide_axis(position = "bottom")) +
  scale_x_continuous(breaks = vbr_df$Station, labels = vbr_df$Station)
vbr_plot

ggsave(vbr_plot, path = "./results/figures/", filename = "VBR_PE_barplot.svg", dpi = 800, height = 3, width = 4.5)



ratio_long <- vp_pe_df %>%
  select(Cruise, Station, VBR, HNALNA, V1V2, V1V3, V2V3) %>%
  pivot_longer(cols = c(VBR, HNALNA, V1V2, V1V3, V2V3),
               names_to = "Ratio",
               values_to = "Count") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),  
         Population = factor(Ratio, levels = c("VBR", "HNALNA", "V1V2", "V1V3", "V2V3"))) 


ratio_grouped_plot <- ggplot(ratio_long, aes(x = Station, y = Count, fill = Ratio)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_hline(yintercept = 1) +
  labs(x = "Station",
       y = "Ratio",
       fill = "Ratio") +
  scale_fill_manual(values = c("VBR" = "#17153B", "HNALNA" = "#1A3636", "V1V2" = "#40534C", "V1V3" = "#677D6A", "V2V3" = "#D6BD98")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside") +
  guides(x = guide_axis(position = "bottom")) +  
  scale_x_continuous(breaks = ratio_long$Station, labels = ratio_long$Station)
ratio_grouped_plot

ggsave(ratio_grouped_plot, path = "./results/figures/", filename = "abundance_ratio_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)



ratio1_long <- vp_pe_df %>%
  select(Cruise, Station, HNALNA, V1V2) %>%
  pivot_longer(cols = c(HNALNA, V1V2),
               names_to = "Ratio",
               values_to = "Count") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),  
         Population = factor(Ratio, levels = c("HNALNA", "V1V2", "V1V3", "V2V3"))) 


hnalna_v1v2_ratio_plot <- ggplot(ratio1_long, aes(x = Station, y = Count, fill = Ratio)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_hline(yintercept = 1) +
  labs(x = "Station",
       y = "Ratio",
       fill = "Ratio") +
  scale_fill_manual(values = c("HNALNA" = "#1A3636", "V1V2" = "#40534C", "V1V3" = "#677D6A", "V2V3" = "#D6BD98")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(Ratio ~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside") +
  guides(x = guide_axis(position = "bottom")) +  
  scale_x_continuous(breaks = ratio_long$Station, labels = ratio_long$Station)
hnalna_v1v2_ratio_plot

ggsave(hnalna_v1v2_ratio_plot, path = "./results/figures/", filename = "hnalna_v1v2_ratio_PE_barplot.svg", dpi = 800, height = 3, width = 4.5)



ratio2_long <- vp_pe_df %>%
  select(Cruise, Station, V1V3, V2V3) %>%
  pivot_longer(cols = c(V1V3, V2V3),
               names_to = "Ratio",
               values_to = "Count") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),  
         Population = factor(Ratio, levels = c("HNALNA", "V1V2", "V1V3", "V2V3"))) 


v1v3_v2v3_ratio_plot <- ggplot(ratio2_long, aes(x = Station, y = Count, fill = Ratio)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_hline(yintercept = 1) +
  labs(x = "Station",
       y = "Ratio",
       fill = "Ratio") +
  scale_fill_manual(values = c("HNALNA" = "#1A3636", "V1V2" = "#40534C", "V1V3" = "#677D6A", "V2V3" = "#D6BD98")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(Ratio ~ Cruise, ncol = 2, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(),
        strip.placement = "outside") +
  guides(x = guide_axis(position = "bottom")) +  
  scale_x_continuous(breaks = ratio_long$Station, labels = ratio_long$Station)

v1v3_v2v3_ratio_plot

ggsave(v1v3_v2v3_ratio_plot, path = "./results/figures/", filename = "v1v3_v2v3_ratio_PE_barplot.svg", dpi = 800, height = 3, width = 4.5)


#### 4.0 Nutrients ####

nutrients <- c("Nitrate", "Nitrite", "Phosphate", "Silicate")
nutrient_colors <- c("Nitrate" = "#1A3636", "Nitrite" = "#40534C", "Phosphate" = "#677D6A", "Silicate" = "#D6BD98")

for (nuts in nutrients) {
  
  nuts_df <- vp_pe_df %>%
    select(Cruise, Station, all_of(nuts)) %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
           Nutrient = nuts)

  nuts_plot <- ggplot(nuts_df, aes(x = Station, y = !!sym(nuts), fill = Nutrient)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Concentration (µM)",
         fill = "Nutrient",
         title = nuts) +
    scale_fill_manual(values = nutrient_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside",
          legend.position = "none") +
    guides(x = guide_axis(position = "bottom")) +
    scale_x_continuous(breaks = nuts_df$Station, labels = nuts_df$Station)
 
  ggsave(nuts_plot, path = "./results/figures/", filename = paste0(nuts, "_PE_barplot.svg"), dpi = 800, height = 3, width = 4.5)
}


#### 5.0 Temperature & Salinity ####


ts_colors <- c("Temperature" = "#973131", "Salinity" = "#E0A75E")

  ts_df <- vp_pe_df %>%
    select(Cruise, Station, Temperature, Salinity) %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")))
  
  temperature_plot <- ggplot(ts_df, aes(x = Station, y = Temperature)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#973131") +
    labs(x = "Station",
         y = "Temperature (°C)",
         title = "Temperature") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +
    scale_x_continuous(breaks = ts_df$Station, labels = ts_df$Station)
  
  ggsave(temperature_plot, path = "./results/figures/", filename = "Temperature_PE_barplot.svg", dpi = 800, height = 3, width = 4.5)


  salinity_plot <- ggplot(ts_df, aes(x = Station, y = Salinity)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#E0A75E") +
    labs(x = "Station",
         y = "Salinity (psu)",
         title = "Salinity") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +
    scale_x_continuous(breaks = ts_df$Station, labels = ts_df$Station)
  
  ggsave(salinity_plot, path = "./results/figures/", filename = "Salinity_PE_barplot.svg", dpi = 800, height = 3, width = 4.5)
  

#### 6.0 Viral Production Rates #####
  # 
  # vp_rate_df <- vp_pe_df %>%
  #   select(Cruise, Station, Lytic_VPCL_c_VP, Lysogenic_VPCL_c_VP) %>%
  #   pivot_longer(cols = c(Lytic_VPCL_c_VP, Lysogenic_VPCL_c_VP),
  #                names_to = "Type",
  #                values_to = "VP") %>%
  #   mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
  #          Type = recode(Type, "Lytic_VPCL_c_VP" = "Lytic", "Lysogenic_VPCL_c_VP" = "Lysogenic"),
  #          Type = factor(Type, levels = c("Lytic", "Lysogenic")),
  #          VP = VP/ 1e5)
  # 
  
  vp_rate_df <- vp_pe_df %>%
    select(Cruise, Station, VP_Lytic , VP_Lysogenic ) %>%
    pivot_longer(cols = c(VP_Lytic, VP_Lysogenic ),
                 names_to = "Type",
                 values_to = "VP") %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
           Type = recode(Type, "VP_Lytic" = "Lytic", "VP_Lysogenic" = "Lysogenic"),
           Type = factor(Type, levels = c("Lytic", "Lysogenic")),
           VP = VP/ 1e4)
  
  vp_rate_plot<- ggplot(vp_rate_df, aes(x = Station, y = VP, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Viral production rate (x 10^4 VLPs mL-1 h-1)",
         fill = "Type") +
      scale_fill_manual(values = c("Lytic" = "#0c1844", "Lysogenic" = "#850000")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = vp_rate_df$Station, labels = vp_rate_df$Station)
  
  vp_rate_plot
  
  ggsave(vp_rate_plot, path = "./results/figures/", filename = "Viral_Production_Rate_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)
  
#### 7.0 Percent lytically infected cells & Lysogeny ####
  
  # Recaluclating using viral production rate
  burst_size <- 25
  vp_pe_df <- vp_pe_df %>%
    mutate(percent_lytic = (VP_Lytic / (Total_Bacteria * burst_size)) * 24 * 100,
           percent_lysogeny = (VP_Lysogenic / (Total_Bacteria * burst_size)) *24 * 100)
  
percent_infection_df <- vp_pe_df %>%
  select(Cruise, Station, percent_lytic, percent_lysogeny) %>%
  pivot_longer(cols = c(percent_lytic, percent_lysogeny),
               names_to = "Type",
               values_to = "Percent_Cells") %>%
  mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
         Type = recode(Type, "percent_lytic" = "Lytically infected", "percent_lysogeny" = "Lysogen"),
         Type = factor(Type, levels = c("Lytically infected", "Lysogen"))#,           VP = VP/ 1e5
         )


  # percent_infection_df <- vp_pe_df %>%
  #   select(Cruise, Station, Lytic_VPCL_Percent_Cells_BS_25, Lysogenic_VPCL_Percent_Cells_BS_25) %>%
  #   pivot_longer(cols = c(Lytic_VPCL_Percent_Cells_BS_25, Lysogenic_VPCL_Percent_Cells_BS_25),
  #                names_to = "Type",
  #                values_to = "Percent_Cells") %>%
  #   mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
  #          Type = recode(Type, "Lytic_VPCL_Percent_Cells_BS_25" = "Lytically infected", "Lysogenic_VPCL_Percent_Cells_BS_25" = "Lysogen"),
  #          Type = factor(Type, levels = c("Lytically infected", "Lysogen"))#,           VP = VP/ 1e5
  #          )
  
  
  percent_infection_plot<- ggplot(percent_infection_df, aes(x = Station, y = Percent_Cells, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Percent bacterial cells (%)",
         fill = "Type") +
    scale_fill_manual(values = c("Lytically infected" = "#0c1844", "Lysogen" = "#850000")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = percent_infection_df$Station, labels = percent_infection_df$Station)
  
  percent_infection_plot
  
  ggsave(percent_infection_plot, path = "./results/figures/", filename = "Percent_lytically_infected_lysogen_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)
  
  
  #### 8.0 Daily bacterial mortality ####
  
  
  
  bacterial_mortality_plot<- ggplot(vp_pe_df, aes(x = Station, y = Lytic_VPCL_Percent_Bacteria_Loss_BS_25)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#0c1844") +
    labs(x = "Station",
         y = "Daily bacterial mortality (%)",
         title = "Daily viral-mediated bacterial mortaity") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = vp_pe_df$Station, labels = vp_pe_df$Station)
  
  bacterial_mortality_plot
  
  ggsave(bacterial_mortality_plot, path = "./results/figures/", filename = "Daily_bacterial_mortality_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 4.5)

  
  #### 9.0 Viral turnover time ####
  
  
  
  viral_turnover_rate_plot<- ggplot(vp_pe_df, aes(x = Station, y = Lytic_LM_V_Turnover_Time)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#0c1844") +
    labs(x = "Station",
         y = "Viral turnover rate (h-1)",
         title = "Viral turnover rate") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = vp_pe_df$Station, labels = vp_pe_df$Station)
  
  viral_turnover_rate_plot
  
  ggsave(viral_turnover_rate_plot, path = "./results/figures/", filename = "Viral_Turnover_Rate_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 4.5)
  
  #### 10.0 Bacterial net growth rate ####
  
  net_growth_rate_df <- vp_pe_df %>%
    select(Cruise, Station, contains("Net_Growth_Rate_")) %>%
    pivot_longer(cols =  contains("Net_Growth_Rate_c_"),
                 names_to = "Population",
                 values_to = "Growth_Rate",
                 names_prefix = "Net_Growth_Rate_c_") %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
          Population = recode(Population, "Bacteria" = "Total_Bacteria"),
           Population = factor(Population, levels = c("Total_Bacteria", "HNA", "LNA")))
  
  
  net_growth_rate_plot<- ggplot(net_growth_rate_df, aes(x = Station, y = Growth_Rate, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Bacterial net growth rate (cells mL-1 h-1)",
         fill = "Population") +
    scale_fill_manual(values = c("Total_Bacteria" = "#850000", "HNA" = "#B80000", "LNA" = "#FF6969")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = net_growth_rate_df$Station, labels = net_growth_rate_df$Station)
  
  net_growth_rate_plot
  
  ggsave(net_growth_rate_plot, path = "./results/figures/", filename = "Bacterial_net_growth_rate_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)
  
  
  
  #### 11.0 Viral decay rate and percent ####
  
  viral_decay_rate_df <- vp_pe_df %>%
    select(Cruise, Station, contains("decay_rate")) %>%
    pivot_longer(cols =  contains("decay_rate_c_"),
                 names_to = "Population",
                 values_to = "Decay_Rate",
                 names_prefix = "decay_rate_c_") %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
           Population = recode(Population, "Viruses" = "Total_Viruses"),
           Population = factor(Population, levels = c("Total_Viruses", "V1", "V2", "V3")),
           Decay_Rate = Decay_Rate/ 1e5)
  
  
  viral_decay_rate_plot<- ggplot(viral_decay_rate_df, aes(x = Station, y = Decay_Rate, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Viral decay rate (x 10 ^5 VLPs mL-1 h-1)",
         fill = "Population") +
    scale_fill_manual(values = c("Total_Viruses" = "#0c1844", "V3" = "#4d778b", "V2" = "#6baed6", "V1" = "#474F7A")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = viral_decay_rate_df$Station, labels = viral_decay_rate_df$Station)
  
  viral_decay_rate_plot
  
  ggsave(viral_decay_rate_plot, path = "./results/figures/", filename = "Viral_decay_rate_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)
  
  
  viral_percent_decay_df <- vp_pe_df %>%
    select(Cruise, Station, contains("percent_decay")) %>%
    pivot_longer(cols =  contains("percent_decay_day_c_"),
                 names_to = "Population",
                 values_to = "Percent_decay",
                 names_prefix = "percent_decay_day_c_") %>%
    mutate(Cruise = factor(Cruise, levels = c("PE477", "PE486")),
           Population = recode(Population, "Viruses" = "Total_Viruses"),
           Population = factor(Population, levels = c("Total_Viruses", "V1", "V2", "V3")))
  
  
  viral_percent_decay_plot<- ggplot(viral_percent_decay_df, aes(x = Station, y = Percent_decay, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
    labs(x = "Station",
         y = "Viral percent decay over day (%)",
         fill = "Population") +
    scale_fill_manual(values = c("Total_Viruses" = "#0c1844", "V3" = "#4d778b", "V2" = "#6baed6", "V1" = "#474F7A")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    facet_wrap(~ Cruise, ncol = 2, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = viral_percent_decay_df$Station, labels = viral_percent_decay_df$Station)
  
  viral_percent_decay_plot
  
  ggsave(viral_percent_decay_plot, path = "./results/figures/", filename = "Viral_percent_decay_day_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 9)
  
  
  
  #### 12.0 Bacterial grazing by heteronanoflagellates ####
  
  
  flb_loss_df <- vp_pe_df %>%
    select(Cruise, Station, contains("flb")) %>%
    mutate(flb_loss_mean_rate = flb_loss_mean_rate/1e3)
  
  
  flb_loss_rate_plot<- ggplot(flb_loss_df, aes(x = Station, y = flb_loss_mean_rate)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#850000") +
    labs(x = "Station",
         y = "FLB loss rate (x 10^3 cells mL-1 h-1)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = flb_loss_df$Station, labels = flb_loss_df$Station)
  
  flb_loss_rate_plot
  
  ggsave(flb_loss_rate_plot, path = "./results/figures/", filename = "FLB_loss_rate_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 4.5)
  
  
  
  flb_loss_percent_decay_plot<- ggplot(flb_loss_df, aes(x = Station, y = percent_flb_grazed_day)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black', fill = "#850000") +
    labs(x = "Station",
         y = "FLB percent decay over day (%)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(),
          strip.placement = "outside") +
    guides(x = guide_axis(position = "bottom")) +  
    scale_x_continuous(breaks = flb_loss_df$Station, labels = flb_loss_df$Station)
  
  flb_loss_percent_decay_plot
  
  ggsave(flb_loss_percent_decay_plot, path = "./results/figures/", filename = "FLB_loss_percent_decay_day_PE_grouped_barplot.svg", dpi = 800, height = 3, width = 4.5)
  
  
  
  