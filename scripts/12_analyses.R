# Setting up ####

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv")

# 2.0 Assigning vectors for parameters ####

physicochemical_params <- c("Temperature", "Salinity", #"Density", "Conductivity", 
                            "Turbidity", "Nitrate", "Phosphate", "Silicate")

biological_params <- c("Oxygen", "Fluorescence",
                "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
                "Total_Viruses", "V1", "V2", "V3",
                "VBR")

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

# 3.0 Correlation matrix for phsyicochemical parameters only####
# To visualize correlations and remove highly correlated variables 

correlation_pc_df <- pe_df %>%
  dplyr::select(all_of(physicochemical_params))

# I chose "pairwise.complete.obs" to remove NAs on a pairwise basis for each variable pairs.

cor_matrix_pc <- cor(correlation_pc_df, use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix_pc <- reshape2::melt(cor_matrix_pc)

# Plot the heatmap
cor_plot_pc<- ggplot(data = melted_cor_matrix_pc, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical & chemical at 3 depths",
       x = NULL,
       y = NULL)
cor_plot_pc

ggsave(plot = cor_plot_pc, path = "./figures/", filename = "correlation_plot_physicochemical_PE_Cruises.svg", dpi = 800, width = 5, height = 5)


# Finding highly correlates variable pairs
highly_correlated_pairs_pc <- melted_cor_matrix_pc %>%
  dplyr::filter(abs(value) > 0.95 & Var1 != Var2) %>%
  arrange(desc(abs(value)))   

print(highly_correlated_pairs_pc)

# Temperature, conductivity, and density are highly correlated at > 0.98. keeping only temperature.
# Turbidity is also highly correlated with temperature, density, and conductivity at > 0.9, but will keep it for now.

correlation_df_low_corr_pc <- correlation_pc_df %>%
  dplyr::select(-c("Density", "Conductivity"))

cor_matrix_low_corr_pc <- cor(correlation_df_low_corr_pc, use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix_low_corr_pc <- reshape2::melt(cor_matrix_low_corr_pc)

# Plot the heatmap
cor_plot_low_corr_pc<- ggplot(data = melted_cor_matrix_low_corr_pc, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical & chemical at 3 depths (low corr)",
       x = NULL,
       y = NULL)
cor_plot_low_corr_pc

ggsave(plot = cor_plot_low_corr_pc, path = "./figures/", filename = "correlation_plot_physicochemical_LOW_CORR_PE_Cruises.svg", dpi = 800, width = 5, height = 5)

correlation_df_low_corr_with_sample_labels_pc <- cbind(pe_df %>% dplyr::select(all_of(sample_params)),
                                                                            correlation_df_low_corr_pc)
  

# 4.0 PCA - low correlated - physicochemical parameters ####

pca_pc_df <- correlation_df_low_corr_with_sample_labels_pc %>%
  na.omit()

pca_pc_labels <- pca_pc_df %>%
  dplyr::select(all_of(sample_params))

pca_pc_numeric_df <- pca_pc_df %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_results <- PCA(pca_pc_numeric_df, scale.unit = T, graph = F)
pca_pc_scores <- as.data.frame(pca_pc_results$ind$coord)
pca_pc_df_final <- cbind(pca_pc_labels, pca_pc_scores)

# Visualizing PCA results
explained_variance_pc <- pca_pc_results$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_loadings <- as.data.frame(pca_pc_results$var$coord)
pca_pc_loadings$Variable <- row.names(pca_pc_loadings)
pca_pc_loadings <- pca_pc_loadings %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc<- ggplot(pca_pc_df_final , aes(x = Dim.1, y = Dim.2, 
                                    color = as.factor(Station_Number))) +  
          stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
          geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
          geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
          scale_color_manual(values = custom_color_palette_stations,
                             guide = guide_legend(ncol = 2)) +
          scale_fill_manual(values = custom_color_palette_cruise) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
          geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
          ggplot2::geom_segment(data = pca_pc_loadings,
                                aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                                arrow = arrow(length = unit(0.2, "cm")),
                                color = "black", linewidth = 1) +  
          geom_text(data = pca_pc_loadings,
                          aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                              y = adjust_label_position(Dim.2, 2, 0.1), 
                              label = Variable),
                          size = 5, color = "black") +  
          scale_shape_manual(values = custom_shape_palette_depths) +  
          labs(
            #title = "PCA Biplot - Physicochemical parameters - 3 depths",
               x = "PC1 (41.35 %)",
               y = "PC2 (38.01 %)",
               color = "Station Number",
               shape = "Depth",
               fill = "Season") +
          theme_test(base_size = 15) +
          theme(legend.position = "right")
pca_pc

ggsave(plot = pca_pc, filename = "./figures/PCA_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)



# 4.1 PCA - PE477 - low correlated - physicochemical parameters ####

pca_pc_1_df <- correlation_df_low_corr_with_sample_labels %>%
  na.omit() %>%
  dplyr::filter(Location == "PE477")

pca_pc_1_labels <- pca_pc_1_df %>%
  dplyr::select(all_of(sample_params))

pca_pc_1_numeric_df <- pca_pc_1_df %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_1_results <- PCA(pca_pc_1_numeric_df, scale.unit = T, graph = F)
pca_pc_1_scores <- as.data.frame(pca_pc_1_results$ind$coord)
pca_pc_1_df_final <- cbind(pca_pc_1_labels, pca_pc_1_scores)

# Visualizing PCA results
explained_variance_pc_1 <- pca_pc_1_results$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_1, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_1), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - PE477",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_1_loadings <- as.data.frame(pca_pc_1_results$var$coord)
pca_pc_1_loadings$Variable <- row.names(pca_pc_1_loadings)
pca_pc_1_loadings <- pca_pc_1_loadings %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_1<- ggplot(pca_pc_1_df_final , aes(x = Dim.1, y = Dim.2, 
                                      color = as.factor(Station_Number))) +  
  #stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_1_loadings,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_1_loadings,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot -PE477 - Physicochemical parameters - 3 depths",
    x = "PC1 (53.43 %)",
    y = "PC2 (22.53 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 20) +
  theme(legend.position = "right")
pca_pc_1

ggsave(plot = pca_pc_1, filename = "./figures/PCA_PE477_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)


# 4.2 PCA - PE486 - low correlated - physicochemical parameters ####

pca_pc_2_df <- correlation_df_low_corr_with_sample_labels %>%
  na.omit() %>%
  dplyr::filter(Location == "PE486")

pca_pc_2_labels <- pca_pc_2_df %>%
  dplyr::select(all_of(sample_params))

pca_pc_2_numeric_df <- pca_pc_2_df %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_2_results <- PCA(pca_pc_2_numeric_df, scale.unit = T, graph = F)
pca_pc_2_scores <- as.data.frame(pca_pc_2_results$ind$coord)
pca_pc_2_df_final <- cbind(pca_pc_2_labels, pca_pc_2_scores)

# Visualizing PCA results
explained_variance_pc_2 <- pca_pc_2_results$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_2, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_2), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - PE486",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_2_loadings <- as.data.frame(pca_pc_2_results$var$coord)
pca_pc_2_loadings$Variable <- row.names(pca_pc_2_loadings)
pca_pc_2_loadings <- pca_pc_2_loadings %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_2<- ggplot(pca_pc_2_df_final , aes(x = Dim.1, y = Dim.2, 
                                          color = as.factor(Station_Number))) +  
  #stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_2_loadings,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_2_loadings,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE486 - Physicochemical parameters - 3 depths",
    x = "PC1 (52.52 %)",
    y = "PC2 (31.39 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 20) +
  theme(legend.position = "right")
pca_pc_2

ggsave(plot = pca_pc_2, filename = "./figures/PCA_PE486_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)


##### 5.0 Box plots and hypothesis testing ####
# 5.1 Physicochemical parameters ####
physicochemical_3depths_df <- pe_df %>%
  dplyr::select(all_of(c(physicochemical_params, sample_params))) %>%
  pivot_longer(cols = all_of(physicochemical_params),
               names_to = "Physicochemical parameters",
               values_to = "Value") %>%
  dplyr::mutate(`Physicochemical parameters` = factor(`Physicochemical parameters`, levels = c("Temperature", "Salinity", "Conductivity", "Density", 
                                                                                             "Turbidity", "Nitrate", "Phosphate", "Silicate")))

# Defining labels using parse-able expressions
pc_variable_labels <- c(
  "Temperature" = "Temperature~(degree*C)",
  "Salinity" = "Salinity~(PSU)",
  "Turbidity" = "Turbidity~(NTU)",
  "Conductivity" = "Conductivity~(mS/cm)",
  "Density" = "Density~(kg/m^3)",  
  "Nitrate" = "Nitrate~(mu*M)",
  "Phosphate" = "Phosphate~(mu*M)",
  "Silicate" = "Silicate~(mu*M)"
)

pc_3d_boxplots <- ggplot(physicochemical_3depths_df,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Physicochemical parameters`, scales = "free_y", nrow = 2, 
             labeller = as_labeller(pc_variable_labels, label_parsed) #, strip.position = "bottom"
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

pc_3d_boxplots

ggsave(plot = pc_3d_boxplots, filename = "./figures/boxplots_physicochemical_3depths.svg", dpi = 800, width = 10, height = 6)

# Perform Kruskal-Wallis test for each physicochemical parameter

kw_results_pc_3d <- physicochemical_3depths_df %>%
  group_by(`Physicochemical parameters`) %>%
  kruskal_test(Value ~ Season) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`Physicochemical parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `Physicochemical parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_pc_3d)

# 5.2 Biological parameters ####
biological_3depths_df <- pe_df %>%
  dplyr::select(all_of(c(biological_params, sample_params))) %>%
  dplyr::mutate(Total_Bacteria = Total_Bacteria/1e+6,
                Total_Viruses = Total_Viruses/1e+6,
                HNA = HNA/1e+6,
                LNA = LNA/1e+6,
                V1 = V1/1e+6,
                V2 = V2/1e+6,
                V3 = V3/1e+6,
                Cyanobacteria = Cyanobacteria/1e+3
                )%>%
  pivot_longer(cols = all_of(biological_params),
               names_to = "Biological parameters",
               values_to = "Value") %>%
  dplyr::mutate(`Biological parameters` = factor(`Biological parameters`, levels = c("Total_Bacteria", "Total_Viruses", "VBR", 
                                                                                     "HNA", "V1","Oxygen",
                                                                                     "LNA", "V2", "Fluorescence",
                                                                                     "Cyanobacteria", "V3")))

# Defining labels using parse-able expressions
bio_variable_labels <- c(
  "Total_Bacteria" = "Bacteria~(10^6~cells~mL^{-1})",
  "Total_Viruses" = "Viruses~(10^6~VLPs~mL^{-1})",
  "VBR" = "VBR",
  "Oxygen" = "Oxygen~(mu*M)",  
  "Fluorescence" = "Fluorescence~(mu*g~L^{-1})",
  "HNA" = "HNA~(10^6~cells~mL^{-1})",
  "LNA" = "LNA~(10^6~cells~mL^{-1})",
  "Cyanobacteria" = "Cyanobacteria~(10^3~cells~mL^{-1})",
  "V1" = "V1~(10^6~VLPs~mL^{-1})",
  "V2" = "V2~(10^6~VLPs~mL^{-1})",
  "V3" = "V3~(10^6~VLPs~mL^{-1})"
)

bio_3d_boxplots <- ggplot(biological_3depths_df,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Biological parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(bio_variable_labels, label_parsed) #, strip.position = "bottom"
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

bio_3d_boxplots

ggsave(plot = bio_3d_boxplots, filename = "./figures/boxplots_biological_3depths.svg", dpi = 800, width = 9, height = 12)



bio_selected_3d_boxplots <- ggplot(biological_3depths_df %>%
                            dplyr::filter(`Biological parameters` %in% c("Total_Bacteria", "Total_Viruses", "VBR",
                                                                      "Cyanobacteria", "Fluorescence", "Oxygen")),
                          aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Biological parameters`, scales = "free_y", nrow = 2, 
             labeller = as_labeller(bio_variable_labels, label_parsed) #, strip.position = "bottom"
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

bio_selected_3d_boxplots

ggsave(plot = bio_selected_3d_boxplots, filename = "./figures/boxplots_biological_selected_3depths.svg", dpi = 800, width = 9, height = 6)





# Perform Kruskal-Wallis test for each biological parameter

kw_results_bio_3d <- biological_3depths_df %>%
  group_by(`Biological parameters`) %>%
  kruskal_test(Value ~ Season) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`Biological parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `Biological parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_bio_3d)


#### 6.0 Testing at 7m depth - PCA & boxplots####

# PCA with all variables in #

pca_pc_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) 


pca_pc_labels_7m <- pca_pc_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_pc_numeric_df_7m <- pca_pc_df_7m %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

pca_pc_results_7m <- PCA(pca_pc_numeric_df_7m, scale.unit = T, graph = F)
pca_pc_scores_7m <- as.data.frame(pca_pc_results_7m$ind$coord)
pca_pc_df_final_7m <- cbind(pca_pc_labels_7m, pca_pc_scores_7m)

# Visualizing PCA results
explained_variance_pc_7m <- pca_pc_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_7m, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")


# PCA biplot

pca_pc_loadings_7m <- as.data.frame(pca_pc_results_7m$var$coord)
pca_pc_loadings_7m$Variable <- row.names(pca_pc_loadings_7m)
pca_pc_loadings_7m <- pca_pc_loadings_7m %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_7m<- ggplot(pca_pc_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                            color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - Physicochemical parameters - 3 depths",
    x = "PC1 (52.19 %)",
    y = "PC2 (30.87 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_pc_7m

ggsave(plot = pca_pc_7m, filename = "./figures/PCA_physicochemical_7m_correlated_not_removed.svg", dpi = 800, width = 10, height = 7)

# Checking for correlation in 7 m physicochemical data ####

correlation_7m_df <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(physicochemical_params))

# I chose "pairwise.complete.obs" to remove NAs on a pairwise basis for each variable pairs.

cor_matrix_7m <- cor(correlation_7m_df, use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix_7m <- reshape2::melt(cor_matrix_7m)

# Plot the heatmap
cor_plot_7m<- ggplot(data = melted_cor_matrix_7m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical & chemical at 7m",
       x = NULL,
       y = NULL)
cor_plot_7m

ggsave(plot = cor_plot_7m, path = "./figures/", filename = "correlation_plot_7m_physicochemical_PE_Cruises.svg", dpi = 800, width = 5, height = 5)


# Finding highly correlates variable pairs
highly_correlated_pairs_7m <- melted_cor_matrix_7m %>%
  dplyr::filter(abs(value) > 0.95 & Var1 != Var2) %>%
  arrange(desc(abs(value)))   

print(highly_correlated_pairs_7m)

# Temperature, conductivity, and density are highly correlated at > 0.98. keeping only temperature.
# Turbidity is also highly correlated with temperature at > 0.95, removing turbidity too

correlation_df_low_corr_7m <- correlation_7m_df %>%
  dplyr::select(-c("Density", "Conductivity", "Turbidity"))

cor_matrix_low_corr_7m <- cor(correlation_df_low_corr_7m, use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix_low_corr_7m <- reshape2::melt(cor_matrix_low_corr_7m)

# Plot the heatmap
cor_plot_low_corr_7m<- ggplot(data = melted_cor_matrix_low_corr_7m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical & chemical at 7m (low corr)",
       x = NULL,
       y = NULL)
cor_plot_low_corr_7m

ggsave(plot = cor_plot_low_corr_7m, path = "./figures/", filename = "correlation_plot_physicochemical_7m_LOW_CORR_PE_Cruises.svg", dpi = 800, width = 5, height = 5)

correlation_7m_df_low_corr_with_sample_labels <- cbind(pe_df %>%
                                                         dplyr::filter(Depth == 7)%>% 
                                                         dplyr::select(all_of(sample_params)),
                                                    correlation_df_low_corr_7m)




#### 6.1 Physicochemical paramters PCA ####

pca_pc_df_7m <- correlation_7m_df_low_corr_with_sample_labels %>%
  na.omit()

pca_pc_labels_7m <- pca_pc_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_pc_numeric_df_7m <- pca_pc_df_7m %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_results_7m <- PCA(pca_pc_numeric_df_7m, scale.unit = T, graph = F)
pca_pc_scores_7m <- as.data.frame(pca_pc_results_7m$ind$coord)
pca_pc_df_final_7m <- cbind(pca_pc_labels_7m, pca_pc_scores_7m)

# Visualizing PCA results
explained_variance_pc_7m <- pca_pc_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_7m, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_loadings_7m <- as.data.frame(pca_pc_results_7m$var$coord)
pca_pc_loadings_7m$Variable <- row.names(pca_pc_loadings_7m)
pca_pc_loadings_7m <- pca_pc_loadings_7m %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_7m<- ggplot(pca_pc_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                      color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - Physicochemical parameters - 3 depths",
    x = "PC1 (52.19 %)",
    y = "PC2 (30.87 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_pc_7m

ggsave(plot = pca_pc_7m, filename = "./figures/PCA_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)

# 6.2 PCA - PE477 - low correlated - physicochemical parameters - 7m ####

pca_pc_1_df_7m <- correlation_7m_df_low_corr_with_sample_labels %>%
  na.omit() %>%
  dplyr::filter(Location == "PE477")

pca_pc_1_labels_7m <- pca_pc_1_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_pc_1_numeric_df_7m <- pca_pc_1_df_7m %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_1_results_7m <- PCA(pca_pc_1_numeric_df_7m, scale.unit = T, graph = F)
pca_pc_1_scores_7m <- as.data.frame(pca_pc_1_results_7m$ind$coord)
pca_pc_1_df_final_7m <- cbind(pca_pc_1_labels_7m, pca_pc_1_scores_7m)

# Visualizing PCA results
explained_variance_pc_1_7m <- pca_pc_1_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_1_7m, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_1_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - PE477 - 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_1_loadings_7m <- as.data.frame(pca_pc_1_results_7m$var$coord)
pca_pc_1_loadings_7m$Variable <- row.names(pca_pc_1_loadings_7m)
pca_pc_1_loadings_7m <- pca_pc_1_loadings_7m %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_1_7m<- ggplot(pca_pc_1_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                          color = as.factor(Station_Number))) +  
  #stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_1_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_1_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot -PE477 - Physicochemical parameters - 7m",
    x = "PC1 (51.18 %)",
    y = "PC2 (35.60 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 20) +
  theme(legend.position = "right")
pca_pc_1_7m

ggsave(plot = pca_pc_1_7m, filename = "./figures/PCA_PE477_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)


# 6.3 PCA - PE486 - low correlated - physicochemical parameters - 7m####

pca_pc_2_df_7m <- correlation_7m_df_low_corr_with_sample_labels %>%
  na.omit() %>%
  dplyr::filter(Location == "PE486")

pca_pc_2_labels_7m <- pca_pc_2_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_pc_2_numeric_df_7m <- pca_pc_2_df_7m %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_2_results_7m <- PCA(pca_pc_2_numeric_df_7m, scale.unit = T, graph = F)
pca_pc_2_scores_7m <- as.data.frame(pca_pc_2_results_7m$ind$coord)
pca_pc_2_df_final_7m <- cbind(pca_pc_2_labels_7m, pca_pc_2_scores_7m)

# Visualizing PCA results
explained_variance_pc_2_7m <- pca_pc_2_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_2_7m, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_2_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance - PE486",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_2_loadings_7m <- as.data.frame(pca_pc_2_results_7m$var$coord)
pca_pc_2_loadings_7m$Variable <- row.names(pca_pc_2_loadings_7m)
pca_pc_2_loadings_7m <- pca_pc_2_loadings_7m %>%
  mutate(Label_X = adjust_label_position(Dim.1, 2, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 2, 0.2))


pca_pc_2_7m <- ggplot(pca_pc_2_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                          color = as.factor(Station_Number))) +  
  #stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_2_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_2_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE486 - Physicochemical parameters - 7m",
    x = "PC1 (65.85 %)",
    y = "PC2 (24.81 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 20) +
  theme(legend.position = "right")
pca_pc_2_7m

ggsave(plot = pca_pc_2_7m, filename = "./figures/PCA_PE486_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)

# 6.4 Physicochemical parameters - 7m ####
physicochemical_7m_df <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(physicochemical_params, sample_params))) %>%
  pivot_longer(cols = all_of(physicochemical_params),
               names_to = "Physicochemical parameters",
               values_to = "Value") %>%
  dplyr::mutate(`Physicochemical parameters` = factor(`Physicochemical parameters`, levels = c("Temperature", "Salinity", "Conductivity", "Density", 
                                                                                               "Turbidity", "Nitrate", "Phosphate", "Silicate")))

# Defining labels using parse-able expressions
pc_variable_labels <- c(
  "Temperature" = "Temperature~(degree*C)",
  "Salinity" = "Salinity~(PSU)",
  "Turbidity" = "Turbidity~(NTU)",
  "Conductivity" = "Conductivity~(mS/cm)",
  "Density" = "Density~(kg/m^3)",  
  "Nitrate" = "Nitrate~(mu*M)",
  "Phosphate" = "Phosphate~(mu*M)",
  "Silicate" = "Silicate~(mu*M)"
)

pc_7m_boxplots <- ggplot(physicochemical_7m_df,
                         aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Physicochemical parameters`, scales = "free_y", nrow = 2, 
             labeller = as_labeller(pc_variable_labels, label_parsed) #, strip.position = "bottom"
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

pc_7m_boxplots

ggsave(plot = pc_7m_boxplots, filename = "./figures/boxplots_physicochemical_7m.svg", dpi = 800, width = 10, height = 6)

# Perform Kruskal-Wallis test for each physicochemical parameter

kw_results_pc_7m <- physicochemical_7m_df %>%
  group_by(`Physicochemical parameters`) %>%
  kruskal_test(Value ~ Season) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`Physicochemical parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `Physicochemical parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_pc_7m)

# 6.5 Biological parameters - 7m ####
biological_7m_df <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(biological_params, sample_params))) %>%
  dplyr::mutate(Total_Bacteria = Total_Bacteria/1e+6,
                Total_Viruses = Total_Viruses/1e+6,
                HNA = HNA/1e+6,
                LNA = LNA/1e+6,
                V1 = V1/1e+6,
                V2 = V2/1e+6,
                V3 = V3/1e+6,
                Cyanobacteria = Cyanobacteria/1e+3
  )%>%
  pivot_longer(cols = all_of(biological_params),
               names_to = "Biological parameters",
               values_to = "Value") %>%
  dplyr::mutate(`Biological parameters` = factor(`Biological parameters`, levels = c("Total_Bacteria", "Total_Viruses", "VBR", 
                                                                                     "HNA", "V1","Oxygen",
                                                                                     "LNA", "V2", "Fluorescence",
                                                                                     "Cyanobacteria", "V3")))

# Defining labels using parse-able expressions
bio_variable_labels <- c(
  "Total_Bacteria" = "Bacteria~(10^6~cells~mL^{-1})",
  "Total_Viruses" = "Viruses~(10^6~VLPs~mL^{-1})",
  "VBR" = "VBR",
  "Oxygen" = "Oxygen~(mu*M)",  
  "Fluorescence" = "Fluorescence~(mu*g~L^{-1})",
  "HNA" = "HNA~(10^6~cells~mL^{-1})",
  "LNA" = "LNA~(10^6~cells~mL^{-1})",
  "Cyanobacteria" = "Cyanobacteria~(10^3~cells~mL^{-1})",
  "V1" = "V1~(10^6~VLPs~mL^{-1})",
  "V2" = "V2~(10^6~VLPs~mL^{-1})",
  "V3" = "V3~(10^6~VLPs~mL^{-1})"
)

bio_7m_boxplots <- ggplot(biological_7m_df,
                          aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Biological parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(bio_variable_labels, label_parsed) #, strip.position = "bottom"
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

bio_7m_boxplots

ggsave(plot = bio_7m_boxplots, filename = "./figures/boxplots_biological_7m.svg", dpi = 800, width = 9, height = 12)



bio_selected_7m_boxplots <- ggplot(biological_7m_df %>%
                                     dplyr::filter(`Biological parameters` %in% c("Total_Bacteria", "Total_Viruses", "VBR",
                                                                                  "Cyanobacteria", "Fluorescence", "Oxygen")),
                                   aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `Biological parameters`, scales = "free_y", nrow = 2, 
             labeller = as_labeller(bio_variable_labels, label_parsed) #, strip.position = "bottom"
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

bio_selected_7m_boxplots

ggsave(plot = bio_selected_7m_boxplots, filename = "./figures/boxplots_biological_selected_7m.svg", dpi = 800, width = 9, height = 6)





# Perform Kruskal-Wallis test for each biological parameter

kw_results_bio_7m <- biological_7m_df %>%
  group_by(`Biological parameters`) %>%
  kruskal_test(Value ~ Season) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`Biological parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `Biological parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_bio_7m)


# 6.6 Viral production parameters - 7m ####
vp_7m_df <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(vp_params, sample_params))) %>%
  dplyr::mutate(decay_rate_linear  = decay_rate_linear /1e+3,
                Corrected_VP_Lytic  = Corrected_VP_Lytic/1e+5,
                Corrected_VP_Lysogenic = Corrected_VP_Lysogenic/1e+5
  )%>%
  pivot_longer(cols = all_of(vp_params),
               names_to = "vp parameters",
               values_to = "Value") %>%
  dplyr::mutate(`vp parameters` = factor(`vp parameters`, levels = c( "Corrected_VP_Lytic", "Corrected_VP_Lysogenic", "decay_rate_linear", 
                                                                      "percent_bacterial_loss_day_burst_50","percent_lysogeny_day_burst_50",
                                                                      "percent_decay_day_linear")))

# Defining labels using parse-able expressions
vp_variable_labels <- c(
  "Corrected_VP_Lytic" = "Lytic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "Corrected_VP_Lysogenic" = "Lysogenic~production~rate~(10^5~cells~mL^{-1}~h^{-1})",
  "decay_rate_linear" = "Viral~decay~rate~(10^3~cells~mL^{-1}~h^{-1})",
  "percent_bacterial_loss_day_burst_50" = "Bacteria~lysed~('%'~cells~d^{-1})",
  "percent_lysogeny_day_burst_50" = "Lysogeny~('%'~cells~d^{-1})",
  "percent_decay_day_linear" = "Viral~decay~('%'~VLPs~d^{-1})"
)

vp_7m_boxplots <- ggplot(vp_7m_df,
                          aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +  
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black", color = "yellow") +  
  facet_wrap(~ `vp parameters`, scales = "free_y", ncol = 3, 
             labeller = as_labeller(vp_variable_labels, label_parsed) #, strip.position = "bottom"
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

vp_7m_boxplots

ggsave(plot = vp_7m_boxplots, filename = "./figures/boxplots_vp_7m.svg", dpi = 800, width = 12, height = 6)




# Perform Kruskal-Wallis test for each VP parameter

kw_results_vp_7m <- vp_7m_df %>%
  group_by(`vp parameters`) %>%
  kruskal_test(Value ~ Season) %>%
  mutate(Significance = case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"  # "ns" for non-significant (p >= 0.05)
  )) %>%
  select(`vp parameters`, statistic, df, p, Significance) %>%
  rename(
    Parameter = `vp parameters`,
    Chi_Square = statistic,
    Degrees_of_Freedom = df,
    P_Value = p
  )
print(kw_results_vp_7m)


# 7.0 PCA with physicochemical, biological, and VP variables at 7M. ####

# Correlations at 7 m ####


correlation_df <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(physicochemical_params, biological_params, vp_params)))

# I chose "pairwise.complete.obs" to remove NAs on a pairwise basis for each variable pairs.

cor_matrix <- cor(correlation_df, use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix <- reshape2::melt(cor_matrix)

# Plot the heatmap
cor_plot<- ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                    hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical, chemical, biological and VP at 7m",
       x = NULL,
       y = NULL)
cor_plot

ggsave(plot = cor_plot, path = "./figures/", filename = "correlation_plot_pc_bio_vp_7m_PE_Cruises.svg", dpi = 800, width = 10, height = 10)


# Finding highly correlates variable pairs
highly_correlated_pairs <- melted_cor_matrix %>%
  dplyr::filter(abs(value) > 0.90 & Var1 != Var2) %>%
  arrange(desc(abs(value)))   

print(highly_correlated_pairs)

# Temperature, conductivity, and density are highly correlated at > 0.98. keeping only temperature.
# Turbidity is also highly correlated with temperature at = 0.95, but will keep it for now.
# Total Viruses and V1 are correlated at > 0.99 and V2 and V1/Total Viruses are correlated at > 0.95
# Total Bacteria and HNA are correlated at > 0.93, but LNA doesn't correlate well. 
# I will do 2 PCAs. One with subgroups and one without. 

correlation_df_main_groups <- correlation_df %>%
  dplyr::select(-c("Density", "Conductivity",
                   "HNA", "LNA", "Cyanobacteria",
                   "V1", "V2", "V3",
                   "decay_rate_linear",
                   contains("Corrected")))

cor_matrix_main_groups  <- cor(correlation_df_main_groups , use = "pairwise.complete.obs", method = "pearson")

melted_cor_matrix_main_groups  <- reshape2::melt(cor_matrix_main_groups)

# Plot the heatmap
cor_plot_main_groups <- ggplot(data = melted_cor_matrix_main_groups , aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: phsyical & chemical at 3 depths (low corr)",
       x = NULL,
       y = NULL)
cor_plot_main_groups 

ggsave(plot = cor_plot_main_groups , path = "./figures/", filename = "correlation_plot_physicochemical_main_groups_PE_Cruises.svg", dpi = 800, width = 5, height = 5)

correlation_df_main_groups_with_sample_labels <- cbind(pe_df %>%
                                                         dplyr::filter(Depth == 7) %>% 
                                                         dplyr::select(all_of(sample_params)),
                                                    correlation_df_main_groups )


pca_df_7m_main_groups <- correlation_df_main_groups_with_sample_labels %>%
  na.omit()

pca_labels_7m_main_groups <- pca_df_7m_main_groups %>%
  dplyr::select(all_of(sample_params))

pca_numeric_df_7m_main_groups <- pca_df_7m_main_groups %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for all main parameters at 7m
pca_results_7m_main_groups <- PCA(pca_numeric_df_7m_main_groups, scale.unit = T, graph = F)
pca_scores_7m_main_groups <- as.data.frame(pca_results_7m_main_groups$ind$coord)
pca_df_final_7m_main_groups <- cbind(pca_labels_7m_main_groups, pca_scores_7m_main_groups)

# Visualizing PCA results
explained_variance_7m_main_groups <- pca_results_7m_main_groups$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_7m_main_groups, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance -Main groups - 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")


# PCA biplot

pca_loadings_7m_main_groups <- as.data.frame(pca_results_7m_main_groups$var$coord)
pca_loadings_7m_main_groups$Variable <- row.names(pca_loadings_7m_main_groups)
pca_loadings_7m_main_groups <- pca_loadings_7m_main_groups %>%
  mutate(Label_X = adjust_label_position(Dim.1, 3, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 3, 0.2))


pca_7m_main_groups <- ggplot(pca_df_final_7m_main_groups , aes(x = Dim.1, y = Dim.2, 
                                            color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  

  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_loadings_7m_main_groups,
                        aes(x = 0, y = 0, xend = Dim.1 * 3, yend = Dim.2 * 3),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_loadings_7m_main_groups,
            aes(x = adjust_label_position(Dim.1, 3, 0.65), 
                y = adjust_label_position(Dim.2, 3, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - All main parameters - 7m",
    x = "PC1 (34.06 %)",
    y = "PC2 (19.06 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_7m_main_groups

ggsave(plot = pca_7m_main_groups, filename = "./figures/PCA_all_main_parameters_7m.svg", dpi = 800, width = 10, height = 7)

# With subgroups #####
correlation_df_sub_groups <- correlation_df %>%
  dplyr::select(-c("Density", "Conductivity",
                   #"HNA", "LNA", "Cyanobacteria",
                   #"V1", "V2", "V3",
                   "decay_rate_linear",
                   contains("Corrected")))

correlation_df_sub_groups_with_sample_labels <- cbind(pe_df %>%
                                                         dplyr::filter(Depth == 7) %>% 
                                                         dplyr::select(all_of(sample_params)),
                                                       correlation_df_sub_groups )


pca_df_7m_sub_groups <- correlation_df_sub_groups_with_sample_labels %>%
  na.omit()

pca_labels_7m_sub_groups <- pca_df_7m_sub_groups %>%
  dplyr::select(all_of(sample_params))

pca_numeric_df_7m_sub_groups <- pca_df_7m_sub_groups %>%
  dplyr::select(-all_of(sample_params)) %>%
  scale()

# Performing PCA for all sub parameters at 7m
pca_results_7m_sub_groups <- PCA(pca_numeric_df_7m_sub_groups, scale.unit = T, graph = F)
pca_scores_7m_sub_groups <- as.data.frame(pca_results_7m_sub_groups$ind$coord)
pca_df_final_7m_sub_groups <- cbind(pca_labels_7m_sub_groups, pca_scores_7m_sub_groups)

# Visualizing PCA results
explained_variance_7m_sub_groups <- pca_results_7m_sub_groups$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_7m_sub_groups, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: PCA Explained Variance -sub groups - 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")


# PCA biplot

pca_loadings_7m_sub_groups <- as.data.frame(pca_results_7m_sub_groups$var$coord)
pca_loadings_7m_sub_groups$Variable <- row.names(pca_loadings_7m_sub_groups)
pca_loadings_7m_sub_groups <- pca_loadings_7m_sub_groups %>%
  mutate(Label_X = adjust_label_position(Dim.1, 3, 0.2),  
         Label_Y = adjust_label_position(Dim.2, 3, 0.2))


pca_7m_sub_groups <- ggplot(pca_df_final_7m_sub_groups , aes(x = Dim.1, y = Dim.2, 
                                                               color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_loadings_7m_sub_groups,
                        aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_loadings_7m_sub_groups,
            aes(x = adjust_label_position(Dim.1, 4, 0.65), 
                y = adjust_label_position(Dim.2, 4, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - All sub parameters - 7m",
    x = "PC1 (37.51 %)",
    y = "PC2 (17.51 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_7m_sub_groups

ggsave(plot = pca_7m_sub_groups, filename = "./figures/PCA_all_sub_parameters_7m.svg", dpi = 800, width = 10, height = 7)

# Biological PCA



#### 7.0 Linear regression models to explain lytic and lysogenic production rates ####
#### 8.0 Dogger's Bank in PE486 ####
# Group 9, 12.1, and 12.2 as Dogger bank stations and 7, 10, 11, 8 as non