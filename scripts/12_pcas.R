# AIM: To perform PCAs #####

# Setting up ####

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv")

# 2.0 Assigning vectors for parameters ####

physicochemical_params <- c("Temperature", "Salinity", #"Density", "Conductivity", 
                            "Turbidity", "Nitrate", "Phosphate", "Silicate")

biological_params <- c(#"Oxygen", 
                       "Chlorophyll",
                       "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
                       "Total_Viruses", "V1", "V2", "V3",
                       "VBR")

vp_params <- c("decay_rate_linear",
              # "percent_decay_day_linear",
               "Corrected_VP_Lytic",
               "Corrected_VP_Lysogenic"
               # "Corrected_VP_SE_Lytic",
               # "Corrected_VP_SE_Lysogenic",
               # "Corrected_VP_Lytic_per_bacteria",
               # "Corrected_VP_Lysogenic_per_bacteria",
               # "percent_bacterial_loss_day_burst_50",
               # "percent_lysogeny_day_burst_50"
               )

sample_params <- c("Location", "Station_Number", "Depth", 
                   "sample_tag", "Latitude", "Longitude",
                   "Season")

# 3.0 PCA - Physicochemical parameters - 3 depths ####

pca_pc_df_3d <- pe_df %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pca_pc_labels_3d <- pca_pc_df_3d %>%
  dplyr::select(all_of(sample_params))

pca_pc_numeric_df_3d <- pca_pc_df_3d %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pca_pc_results_3d <- PCA(pca_pc_numeric_df_3d, scale.unit = T, graph = F)
pca_pc_scores_3d <- as.data.frame(pca_pc_results_3d$ind$coord)
pca_pc_df_final_3d  <- cbind(pca_pc_labels_3d , pca_pc_scores_3d )

# Visualizing PCA results
explained_variance_pc_3d  <- pca_pc_results_3d$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_3d , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc_3d), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_pc_loadings_3d  <- as.data.frame(pca_pc_results_3d$var$coord)
pca_pc_loadings_3d$Variable <- row.names(pca_pc_loadings_3d )


pca_pc_3d <- ggplot(pca_pc_df_final_3d , aes(x = Dim.1, y = Dim.2, 
                                      color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_pc_loadings_3d,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_pc_loadings_3d,
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
pca_pc_3d

ggsave(plot = pca_pc_3d, filename = "./figures/PCA_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)



# 3.1 PE477 - PCA - Physicochemical parameters - 3 depths ####

pe477_pca_pc_df_3d <- pe_df %>%
  dplyr::filter(Location == "PE477") %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pe477_pca_pc_labels_3d <- pe477_pca_pc_df_3d %>%
  dplyr::select(all_of(sample_params))

pe477_pca_pc_numeric_df_3d <- pe477_pca_pc_df_3d %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pe477_pca_pc_results_3d <- PCA(pe477_pca_pc_numeric_df_3d, scale.unit = T, graph = F)
pe477_pca_pc_scores_3d <- as.data.frame(pe477_pca_pc_results_3d$ind$coord)
pe477_pca_pc_df_final_3d  <- cbind(pe477_pca_pc_labels_3d , pe477_pca_pc_scores_3d )

# Visualizing PCA results
pe477_explained_variance_pc_3d  <- pe477_pca_pc_results_3d$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe477_explained_variance_pc_3d , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe477_explained_variance_pc_3d), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")

# PCA biplot

pe477_pca_pc_loadings_3d  <- as.data.frame(pe477_pca_pc_results_3d$var$coord)
pe477_pca_pc_loadings_3d$Variable <- row.names(pe477_pca_pc_loadings_3d )


pe477_pca_pc_3d <- ggplot(pe477_pca_pc_df_final_3d , aes(x = Dim.1, y = Dim.2, 
                                             color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe477_pca_pc_loadings_3d,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe477_pca_pc_loadings_3d,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE477- Physicochemical parameters - 3 depths",
    x = "PC1 (53.43 %)",
    y = "PC2 (22.53 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe477_pca_pc_3d

ggsave(plot = pe477_pca_pc_3d, filename = "./figures/PE477_PCA_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)


# 3.2 PE486 - PCA - Physicochemical parameters - 3 depths ####

pe486_pca_pc_df_3d <- pe_df %>%
  dplyr::filter(Location == "PE486") %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pe486_pca_pc_labels_3d <- pe486_pca_pc_df_3d %>%
  dplyr::select(all_of(sample_params))

pe486_pca_pc_numeric_df_3d <- pe486_pca_pc_df_3d %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 3 depths
pe486_pca_pc_results_3d <- PCA(pe486_pca_pc_numeric_df_3d, scale.unit = T, graph = F)
pe486_pca_pc_scores_3d <- as.data.frame(pe486_pca_pc_results_3d$ind$coord)
pe486_pca_pc_df_final_3d  <- cbind(pe486_pca_pc_labels_3d , pe486_pca_pc_scores_3d )

# Visualizing PCA results
pe486_explained_variance_pc_3d  <- pe486_pca_pc_results_3d$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe486_explained_variance_pc_3d , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe486_explained_variance_pc_3d), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pe486_pca_pc_loadings_3d  <- as.data.frame(pe486_pca_pc_results_3d$var$coord)
pe486_pca_pc_loadings_3d$Variable <- row.names(pe486_pca_pc_loadings_3d )


pe486_pca_pc_3d <- ggplot(pe486_pca_pc_df_final_3d , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe486_pca_pc_loadings_3d,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe486_pca_pc_loadings_3d,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - pe486- Physicochemical parameters - 3 depths",
    x = "PC1 (52.52 %)",
    y = "PC2 (31.39 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe486_pca_pc_3d

ggsave(plot = pe486_pca_pc_3d, filename = "./figures/PE486_PCA_physicochemical_3depths.svg", dpi = 800, width = 10, height = 7)



# 4.0 PCA - Physicochemical parameters - 7m ####

pca_pc_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pca_pc_labels_7m <- pca_pc_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_pc_numeric_df_7m <- pca_pc_df_7m %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 7m
pca_pc_results_7m <- PCA(pca_pc_numeric_df_7m, scale.unit = T, graph = F)
pca_pc_scores_7m <- as.data.frame(pca_pc_results_7m$ind$coord)
pca_pc_df_final_7m  <- cbind(pca_pc_labels_7m , pca_pc_scores_7m )

# Visualizing PCA results
explained_variance_pc_7m  <- pca_pc_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_pc_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_pc), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 80 % of the variance

# PCA biplot

pca_pc_loadings_7m  <- as.data.frame(pca_pc_results_7m$var$coord)
pca_pc_loadings_7m$Variable <- row.names(pca_pc_loadings_7m )


pca_pc_7m <- ggplot(pca_pc_df_final_7m , aes(x = Dim.1, y = Dim.2, 
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
    #title = "PCA Biplot - Physicochemical parameters - 7m",
    x = "PC1 (48.35 %)",
    y = "PC2 (35.65 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_pc_7m

ggsave(plot = pca_pc_7m, filename = "./figures/PCA_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)



# 4.1 PE477 - PCA - Physicochemical parameters - 7 m ####

pe477_pca_pc_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE477",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pe477_pca_pc_labels_7m <- pe477_pca_pc_df_7m %>%
  dplyr::select(all_of(sample_params))

pe477_pca_pc_numeric_df_7m <- pe477_pca_pc_df_7m %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 7 m
pe477_pca_pc_results_7m <- PCA(pe477_pca_pc_numeric_df_7m, scale.unit = T, graph = F)
pe477_pca_pc_scores_7m <- as.data.frame(pe477_pca_pc_results_7m$ind$coord)
pe477_pca_pc_df_final_7m  <- cbind(pe477_pca_pc_labels_7m , pe477_pca_pc_scores_7m )

# Visualizing PCA results
pe477_explained_variance_pc_7m  <- pe477_pca_pc_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe477_explained_variance_pc_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe477_explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 7 m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")

# PCA biplot

pe477_pca_pc_loadings_7m  <- as.data.frame(pe477_pca_pc_results_7m$var$coord)
pe477_pca_pc_loadings_7m$Variable <- row.names(pe477_pca_pc_loadings_7m )


pe477_pca_pc_7m <- ggplot(pe477_pca_pc_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe477_pca_pc_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe477_pca_pc_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE477- Physicochemical parameters - 3 depths",
    x = "PC1 (49.41 %)",
    y = "PC2 (36.67 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe477_pca_pc_7m

ggsave(plot = pe477_pca_pc_7m, filename = "./figures/PE477_PCA_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)


# 4.2 PE486 - PCA - Physicochemical parameters - 7 m ####

pe486_pca_pc_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE486",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params))) %>%
  na.omit()

pe486_pca_pc_labels_7m <- pe486_pca_pc_df_7m %>%
  dplyr::select(all_of(sample_params))

pe486_pca_pc_numeric_df_7m <- pe486_pca_pc_df_7m %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  scale()

# Performing PCA for physicochemical paramters at 7 m
pe486_pca_pc_results_7m <- PCA(pe486_pca_pc_numeric_df_7m, scale.unit = T, graph = F)
pe486_pca_pc_scores_7m <- as.data.frame(pe486_pca_pc_results_7m$ind$coord)
pe486_pca_pc_df_final_7m  <- cbind(pe486_pca_pc_labels_7m , pe486_pca_pc_scores_7m )

# Visualizing PCA results
pe486_explained_variance_pc_7m  <- pe486_pca_pc_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe486_explained_variance_pc_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe486_explained_variance_pc_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pe486_pca_pc_loadings_7m  <- as.data.frame(pe486_pca_pc_results_7m$var$coord)
pe486_pca_pc_loadings_7m$Variable <- row.names(pe486_pca_pc_loadings_7m )


pe486_pca_pc_7m <- ggplot(pe486_pca_pc_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe486_pca_pc_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe486_pca_pc_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - pe486- Physicochemical parameters - 7 m",
    x = "PC1 (56.38 %)",
    y = "PC2 (27.80 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe486_pca_pc_7m

ggsave(plot = pe486_pca_pc_7m, filename = "./figures/PE486_PCA_physicochemical_7m.svg", dpi = 800, width = 10, height = 7)

# 5.0 PCA - Bioogical parameters - 7m ####

pca_bio_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, biological_params))) %>%
  na.omit()

pca_bio_labels_7m <- pca_bio_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_bio_numeric_df_7m <- pca_bio_df_7m %>%
  dplyr::select(all_of(biological_params)) %>%
  scale()

# Performing PCA for biological paramters at 7m
pca_bio_results_7m <- PCA(pca_bio_numeric_df_7m, scale.unit = T, graph = F)
pca_bio_scores_7m <- as.data.frame(pca_bio_results_7m$ind$coord)
pca_bio_df_final_7m  <- cbind(pca_bio_labels_7m , pca_bio_scores_7m )

# Visualizing PCA results
explained_variance_bio_7m  <- pca_bio_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_bio_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_bio_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: biological parameters at 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 80 % of the variance

# PCA biplot

pca_bio_loadings_7m  <- as.data.frame(pca_bio_results_7m$var$coord)
pca_bio_loadings_7m$Variable <- row.names(pca_bio_loadings_7m )


pca_bio_7m <- ggplot(pca_bio_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                             color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_bio_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_bio_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - biological parameters - 7m",
    x = "PC1 (52.99 %)",
    y = "PC2 (24.23 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_bio_7m

ggsave(plot = pca_bio_7m, filename = "./figures/PCA_biological_7m.svg", dpi = 800, width = 10, height = 7)



# 5.1 PE477 - PCA - biological parameters - 7 m ####

pe477_pca_bio_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE477",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, biological_params))) %>%
  na.omit()

pe477_pca_bio_labels_7m <- pe477_pca_bio_df_7m %>%
  dplyr::select(all_of(sample_params))

pe477_pca_bio_numeric_df_7m <- pe477_pca_bio_df_7m %>%
  dplyr::select(all_of(biological_params)) %>%
  scale()

# Performing PCA for biological paramters at 7 m
pe477_pca_bio_results_7m <- PCA(pe477_pca_bio_numeric_df_7m, scale.unit = T, graph = F)
pe477_pca_bio_scores_7m <- as.data.frame(pe477_pca_bio_results_7m$ind$coord)
pe477_pca_bio_df_final_7m  <- cbind(pe477_pca_bio_labels_7m , pe477_pca_bio_scores_7m )

# Visualizing PCA results
pe477_explained_variance_bio_7m  <- pe477_pca_bio_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe477_explained_variance_bio_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe477_explained_variance_bio_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: biological parameters at 7 m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")

# PCA biplot

pe477_pca_bio_loadings_7m  <- as.data.frame(pe477_pca_bio_results_7m$var$coord)
pe477_pca_bio_loadings_7m$Variable <- row.names(pe477_pca_bio_loadings_7m )


pe477_pca_bio_7m <- ggplot(pe477_pca_bio_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe477_pca_bio_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe477_pca_bio_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE477- biological parameters - 3 depths",
    x = "PC1 (49.87 %)",
    y = "PC2 (34.31 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe477_pca_bio_7m

ggsave(plot = pe477_pca_bio_7m, filename = "./figures/PE477_PCA_biological_7m.svg", dpi = 800, width = 10, height = 7)


# 5.2 PE486 - PCA - biological parameters - 7 m ####

pe486_pca_bio_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE486",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, biological_params))) %>%
  na.omit()

pe486_pca_bio_labels_7m <- pe486_pca_bio_df_7m %>%
  dplyr::select(all_of(sample_params))

pe486_pca_bio_numeric_df_7m <- pe486_pca_bio_df_7m %>%
  dplyr::select(all_of(biological_params)) %>%
  scale()

# Performing PCA for biological paramters at 7 m
pe486_pca_bio_results_7m <- PCA(pe486_pca_bio_numeric_df_7m, scale.unit = T, graph = F)
pe486_pca_bio_scores_7m <- as.data.frame(pe486_pca_bio_results_7m$ind$coord)
pe486_pca_bio_df_final_7m  <- cbind(pe486_pca_bio_labels_7m , pe486_pca_bio_scores_7m )

# Visualizing PCA results
pe486_explained_variance_bio_7m  <- pe486_pca_bio_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe486_explained_variance_bio_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe486_explained_variance_bio_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: biological parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pe486_pca_bio_loadings_7m  <- as.data.frame(pe486_pca_bio_results_7m$var$coord)
pe486_pca_bio_loadings_7m$Variable <- row.names(pe486_pca_bio_loadings_7m )


pe486_pca_bio_7m <- ggplot(pe486_pca_bio_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe486_pca_bio_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe486_pca_bio_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - pe486- biological parameters - 7 m",
    x = "PC1 (65.84 %)",
    y = "PC2 (22.66 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe486_pca_bio_7m

ggsave(plot = pe486_pca_bio_7m, filename = "./figures/PE486_PCA_biological_7m.svg", dpi = 800, width = 10, height = 7)


# 5.3 bIOLOGICAL 3 DEPTHS ####

pca_bio_df_3d <- pe_df %>%
  dplyr::select(all_of(c(sample_params, biological_params))) %>%
  na.omit()

pca_bio_labels_3d <- pca_bio_df_3d %>%
  dplyr::select(all_of(sample_params))

pca_bio_numeric_df_3d <- pca_bio_df_3d %>%
  dplyr::select(all_of(biological_params)) %>%
  scale()

# Performing PCA for biological paramters at 3 depths
pca_bio_results_3d <- PCA(pca_bio_numeric_df_3d, scale.unit = T, graph = F)
pca_bio_scores_3d <- as.data.frame(pca_bio_results_3d$ind$coord)
pca_bio_df_final_3d  <- cbind(pca_bio_labels_3d , pca_bio_scores_3d )

# Visualizing PCA results
explained_variance_bio_3d  <- pca_bio_results_3d$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_bio_3d , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_bio_3d), 1)) +
  theme_test() +
  labs(title = "Scree Plot: biological parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pca_bio_loadings_3d  <- as.data.frame(pca_bio_results_3d$var$coord)
pca_bio_loadings_3d$Variable <- row.names(pca_bio_loadings_3d )


pca_bio_3d <- ggplot(pca_bio_df_final_3d , aes(x = Dim.1, y = Dim.2, 
                                             color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_bio_loadings_3d,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_bio_loadings_3d,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - biological parameters - 3 depths",
    x = "PC1 (51.95 %)",
    y = "PC2 (24.84 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_bio_3d

ggsave(plot = pca_bio_3d, filename = "./figures/PCA_biological_3depths.svg", dpi = 800, width = 10, height = 7)





# 6.0 PCA - physicochemical, biological and VP - 7m ####


pca_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params, biological_params, vp_params))) %>%
  na.omit()

pca_labels_7m <- pca_df_7m %>%
  dplyr::select(all_of(sample_params))

pca_numeric_df_7m <- pca_df_7m %>%
  dplyr::select(all_of(c(physicochemical_params, biological_params, vp_params))) %>%
  scale()

# Performing PCA for physicochemical paramters at 7m
pca_results_7m <- PCA(pca_numeric_df_7m, scale.unit = T, graph = F)
pca_scores_7m <- as.data.frame(pca_results_7m$ind$coord)
pca_df_final_7m  <- cbind(pca_labels_7m , pca_scores_7m )

# Visualizing PCA results
explained_variance_7m  <- pca_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(explained_variance_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(explained_variance_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical, bio, & VP parameters at 7m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")

# PCA biplot

pca_loadings_7m  <- as.data.frame(pca_results_7m$var$coord)
pca_loadings_7m$Variable <- row.names(pca_loadings_7m )


pca_7m <- ggplot(pca_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                             color = as.factor(Station_Number))) +  
  stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pca_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 5, yend = Dim.2 * 5),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pca_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 5, 0.65), 
                y = adjust_label_position(Dim.2, 5, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - Physicochemical parameters - 7m",
    x = "PC1 (47.93 %)",
    y = "PC2 (19.47 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pca_7m

ggsave(plot = pca_7m, filename = "./figures/PCA_all_parameters_7m.svg", dpi = 800, width = 10, height = 7)



# 5.1 PE477 - PCA - Physicochemical parameters - 7 m ####

pe477_pca_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE477",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params, biological_params, vp_params))) %>%
  na.omit()

pe477_pca_labels_7m <- pe477_pca_df_7m %>%
  dplyr::select(all_of(sample_params))

pe477_pca_numeric_df_7m <- pe477_pca_df_7m %>%
  dplyr::select(all_of(c(physicochemical_params, biological_params, vp_params))) %>%
  scale()

# Performing PCA for physicochemical paramters at 7 m
pe477_pca_results_7m <- PCA(pe477_pca_numeric_df_7m, scale.unit = T, graph = F)
pe477_pca_scores_7m <- as.data.frame(pe477_pca_results_7m$ind$coord)
pe477_pca_df_final_7m  <- cbind(pe477_pca_labels_7m , pe477_pca_scores_7m )

# Visualizing PCA results
pe477_explained_variance_7m  <- pe477_pca_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe477_explained_variance_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe477_explained_variance_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 7 m",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")

# PCA biplot

pe477_pca_loadings_7m  <- as.data.frame(pe477_pca_results_7m$var$coord)
pe477_pca_loadings_7m$Variable <- row.names(pe477_pca_loadings_7m )


pe477_pca_7m <- ggplot(pe477_pca_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe477_pca_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe477_pca_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - PE477- Physicochemical parameters - 3 depths",
    x = "PC1 (41.36 %)",
    y = "PC2 (27.97 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe477_pca_7m

ggsave(plot = pe477_pca_7m, filename = "./figures/PE477_PCA_all_parameters_7m.svg", dpi = 800, width = 10, height = 7)


# 5.2 PE486 - PCA - Physicochemical parameters - 7 m ####

pe486_pca_df_7m <- pe_df %>%
  dplyr::filter(Location == "PE486",
                Depth == 7) %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params, biological_params, vp_params))) %>%
  na.omit()

pe486_pca_labels_7m <- pe486_pca_df_7m %>%
  dplyr::select(all_of(sample_params))

pe486_pca_numeric_df_7m <- pe486_pca_df_7m %>%
  dplyr::select(all_of(c(physicochemical_params, biological_params, vp_params))) %>%
  scale()

# Performing PCA for physicochemical paramters at 7 m
pe486_pca_results_7m <- PCA(pe486_pca_numeric_df_7m, scale.unit = T, graph = F)
pe486_pca_scores_7m <- as.data.frame(pe486_pca_results_7m$ind$coord)
pe486_pca_df_final_7m  <- cbind(pe486_pca_labels_7m , pe486_pca_scores_7m )

# Visualizing PCA results
pe486_explained_variance_7m  <- pe486_pca_results_7m$eig %>%
  as.data.frame() %>%
  mutate(PC = row_number()) %>%
  rename(Variance = `percentage of variance`)

# Scree plot
ggplot(pe486_explained_variance_7m , aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "black", color = "black") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 4) +
  scale_x_continuous(breaks = seq(1, nrow(pe486_explained_variance_7m), 1)) +
  theme_test() +
  labs(title = "Scree Plot: Physicochemical parameters at 3 depths",
       x = "Principal Component (PC)",
       y = "Variance Explained (%)")
# Together the first 2 PCs explain 79.36 % of the variance

# PCA biplot

pe486_pca_loadings_7m  <- as.data.frame(pe486_pca_results_7m$var$coord)
pe486_pca_loadings_7m$Variable <- row.names(pe486_pca_loadings_7m )


pe486_pca_7m <- ggplot(pe486_pca_df_final_7m , aes(x = Dim.1, y = Dim.2, 
                                                         color = as.factor(Station_Number))) +  
  # stat_ellipse(aes(group = Season, fill = Season), geom = "polygon", alpha = 0.2, level = 0.95) +  
  geom_point(aes(shape = as.factor(Depth)),size = 4, alpha = 0.8) +  
  geom_text_repel(aes(label = Station_Number), color = "black", size = 5, box.padding = 0.3, max.overlaps = 100) + 
  scale_color_manual(values = custom_color_palette_stations,
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = custom_color_palette_cruise) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  ggplot2::geom_segment(data = pe486_pca_loadings_7m,
                        aes(x = 0, y = 0, xend = Dim.1 * 2, yend = Dim.2 * 2),
                        arrow = arrow(length = unit(0.2, "cm")),
                        color = "black", linewidth = 1) +  
  geom_text(data = pe486_pca_loadings_7m,
            aes(x = adjust_label_position(Dim.1, 2, 0.65), 
                y = adjust_label_position(Dim.2, 2, 0.1), 
                label = Variable),
            size = 5, color = "black") +  
  scale_shape_manual(values = custom_shape_palette_depths) +  
  labs(
    #title = "PCA Biplot - pe486- Physicochemical parameters - 7 m",
    x = "PC1 (55.50 %)",
    y = "PC2 (17.37 %)",
    color = "Station Number",
    shape = "Depth",
    fill = "Season") +
  theme_test(base_size = 15) +
  theme(legend.position = "right")
pe486_pca_7m

ggsave(plot = pe486_pca_7m, filename = "./figures/PE486_PCA_all_parameters_7m.svg", dpi = 800, width = 10, height = 7)

