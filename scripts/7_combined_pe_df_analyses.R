library(tidyverse)
library(ggbiplot)
library(ggsci)
library(ggrepel)
library(stats)
library(colorspace)
library(ggforce)

vp_pe_df <- read.csv("./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_analyzed_pe_cruises_combined.txt", sep = "\t")

pe477_mean_temp <- vp_pe_df %>%
  filter(Cruise == 'PE477') %>%
  select(Temperature) %>%
  summarise(mean = mean(Temperature, na.rm = T))
vp_pe_df$Temperature <- vp_pe_df$Temperature %>%
  replace_na(pe477_mean_temp[1,])

pe477_mean_salinity <- vp_pe_df %>%
  filter(Cruise == 'PE477') %>%
  select(Salinity) %>%
    summarise(mean = mean(Salinity, na.rm = T))
vp_pe_df$Salinity <- vp_pe_df$Salinity %>%
  replace_na(pe477_mean_salinity[1,])


# removing flb data as it is only present for one cruise

vp_pe_df <- vp_pe_df %>%
  select(-contains("flb")) %>%
  select(-contains("LM"), 
         # -contains("BS_10"),
         # -contains("BS_40"),
         -contains("tude"),
         -contains("VP_SE"),
         -contains("DON"),
         -contains("DOP"),
         -contains("DOC"),
         -contains("Generation"),
         -c(Station_Number2, Depth) 
         #-c(B_0, B_OS, V_OS)
         ) %>%
  dplyr::mutate(across(contains("decay"), ~ . * -1))
colnames(vp_pe_df) 


#Writing this out for ODV

write.table(vp_pe_df2, "./PE_Cruises/results/PE_Cruises_viral_production/analyze_pe_VIPCAL_LM/vp_analyzed_pe_cruises_odv_input.txt", sep = "\t",  row.names = F, quote = F)




numeric_columns <- vp_pe_df %>%
  select_if(is.numeric)

which(is.na(numeric_columns))


# cHECKIGN FOR COVARIATES
cor_matrix <- cor(numeric_columns)


# Melt the correlation matrix for ggplot
melted_cor_matrix <- reshape2::melt(cor_matrix)

# Plot the heatmap
cor_plot<- ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Heatmap PE Cruises",
       x = NULL,
       y = NULL)

cor_plot

ggsave(plot = cor_plot, path = "./results/figures/", filename = "correlation_plot_PE_Cruises.svg", dpi = 800, width = 10, height = 10)


cor_matrix <- cor(numeric_columns, use = "complete.obs")

# Find indices of highly correlated variables (e.g., correlation > 0.95)
highly_correlated_indices <- caret::findCorrelation(cor_matrix, cutoff = 0.95, names = TRUE)

# Get the names of the highly correlated pairs
highly_correlated_pairs <- combn(highly_correlated_indices, 2, function(x) if(abs(cor_matrix[x[1], x[2]]) > 0.95) x, simplify = FALSE)
highly_correlated_pairs <- Filter(Negate(is.null), highly_correlated_pairs)

# Display the highly correlated pairs
highly_correlated_pairs_df <- do.call(rbind, lapply(highly_correlated_pairs, function(x) data.frame(Var1 = x[1], Var2 = x[2], Correlation = cor_matrix[x[1], x[2]])))
print(highly_correlated_pairs_df)



highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9)

# Remove highly correlated variables
numeric_columns_filtered <- numeric_columns %>% select(-highly_correlated[-4]) %>%
  select(-Station) %>%
  select(-contains("abs")) %>%
  select(-contains("decay")) %>%
  select(-contains("V1")) %>%
  select(-contains("V3")) %>%
  select(-contains("Growth"))


cor_matrix_filtered <- cor(numeric_columns_filtered)

# Melt the correlation matrix for ggplot
melted_cor_matrix_filtered <- reshape2::melt(cor_matrix_filtered)

# Plot the heatmap
cor_filtered_plot<- ggplot(data = melted_cor_matrix_filtered, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Heatmap of Remaining Variables",
       x = "Variables",
       y = "Variables")

cor_filtered_plot


ggsave(plot = cor_filtered_plot, path = "./results/figures/", filename = "correlation_filtered_plot_PE_Cruises.svg", dpi = 800, width = 10, height = 10)


# Prepare the data by selecting numeric columns
numeric_columns <- vp_pe_df %>%
  select_if(is.numeric) %>%
  select_if(~ var(.) != 0)

# Standardize the data
numeric_columns_scaled <- scale(numeric_columns_filtered)



# Perform PCA ####
pca_result <- prcomp(numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
pe_biplot<- ggbiplot(pca_result, 
         obs.scale = 1, 
         var.scale = 2.5, 
         groups = pca_data$Cruise, # Assuming you want to color points by Cruise
         ellipse = F, 
         circle = TRUE,
         varname.adjust = 1.5) +
  geom_mark_ellipse(aes(fill = pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-10, 10) +
  ylim(-10, 10)

pe_biplot

ggsave(pe_biplot, path = "./results/figures/", filename = "pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)


#### PCA Abiotic plot #####

abiotic_cols <- c("Nitrite", "Nitrate", "Silicate", "Phosphate", "Temperature")

abiotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(all_of(abiotic_cols))


# Standardize the data
abiotic_numeric_columns_scaled <- scale(abiotic_numeric_columns_filtered)



# Perform PCA ####
abiotic_pca_result <- prcomp(abiotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(abiotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
abiotic_pca_data <- as.data.frame(abiotic_pca_result$x)
abiotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
abiotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(abiotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
abiotic_pe_biplot<- ggbiplot(abiotic_pca_result, 
                     obs.scale = 1, 
                     var.scale = 2.5, 
                     groups = abiotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                     ellipse = F, 
                     circle = TRUE,
                     varname.adjust = 1.5) +
  geom_mark_ellipse(aes(fill = abiotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = abiotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-3, 6) +
  ylim(-2.5, 2.5)

abiotic_pe_biplot

ggsave(abiotic_pe_biplot, path = "./results/figures/", filename = "abiotic_pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)


#### PCA Abiotic plot without temp #####

abiotic_cols <- c("Nitrite", "Nitrate", "Silicate", "Phosphate")

abiotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(all_of(abiotic_cols))


# Standardize the data
abiotic_numeric_columns_scaled <- scale(abiotic_numeric_columns_filtered)



# Perform PCA ####
abiotic_pca_result <- prcomp(abiotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(abiotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
abiotic_pca_data <- as.data.frame(abiotic_pca_result$x)
abiotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
abiotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(abiotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
abiotic_pe_biplot<- ggbiplot(abiotic_pca_result, 
                             obs.scale = 1, 
                             var.scale = 2.5, 
                             groups = abiotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                             ellipse = F, 
                             circle = TRUE,
                             varname.adjust = 1.5) +
  #geom_mark_ellipse(aes(fill = abiotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = abiotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  # scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  # scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-3, 6) +
  ylim(-2.5, 2.5)

abiotic_pe_biplot

ggsave(abiotic_pe_biplot, path = "./results/figures/", filename = "abiotic_pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)


#### PCA Biotic plot #####

abiotic_cols_to_remove <- c("Nitrite", "Nitrate", "Silicate", "Phosphate", "Salinity",  "Temperature")

biotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(!all_of(abiotic_cols_to_remove))


# Standardize the data
biotic_numeric_columns_scaled <- scale(biotic_numeric_columns_filtered)



# Perform PCA ####
biotic_pca_result <- prcomp(biotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(biotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
biotic_pca_data <- as.data.frame(biotic_pca_result$x)
biotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
biotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(biotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
biotic_pe_biplot<- ggbiplot(biotic_pca_result, 
                             obs.scale = 1, 
                             var.scale = 2.5, 
                             groups = biotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                             ellipse = F, 
                             circle = F,
                             varname.adjust = 1.5) +
  #geom_mark_ellipse(aes(fill = biotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = biotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  # scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  # scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-10, 10) +
  ylim(-10.5, 9.5)


biotic_pe_biplot

ggsave(biotic_pe_biplot, path = "./results/figures/", filename = "biotic_pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)


###mORE EXCLUSIVE

#### PCA Biotic plot #####
cols_to_remove <- c("Nitrite", "Nitrate", "Silicate", "Phosphate", "Salinity", "Temperature", "LNA", "HNALNA")

biotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(-all_of(cols_to_remove), 
                -matches("Growth"), 
                -matches("decay"),
                -matches("abs"), 
                -matches("V1"), 
                -matches("V2"), 
                -matches("V3"))


# Standardize the data
biotic_numeric_columns_scaled <- scale(biotic_numeric_columns_filtered)



# Perform PCA ####
biotic_pca_result <- prcomp(biotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(biotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
biotic_pca_data <- as.data.frame(biotic_pca_result$x)
biotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
biotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(biotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
biotic_pe_biplot<- ggbiplot(biotic_pca_result, 
                            obs.scale = 1, 
                            var.scale = 2.5, 
                            groups = biotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                            ellipse = F, 
                            circle = F,
                            varname.adjust = 1.5) +
  # geom_mark_ellipse(aes(fill = biotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = biotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  # scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  # scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-6, 6) +
  ylim(-8, 7)


biotic_pe_biplot

ggsave(biotic_pe_biplot, path = "./results/figures/", filename = "biotic_pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)





#### PCA Biotic plot #####
cols_to_remove <- c("Salinity",  "LNA", "HNALNA", "Lytic_VPCL_V_Turnover_Time")

biotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(-all_of(cols_to_remove), 
                -matches("Growth"), 
                -matches("decay"), 
                -matches("V1"), 
                -matches("V2"), 
                -matches("V3"))


# Standardize the data
biotic_numeric_columns_scaled <- scale(biotic_numeric_columns_filtered)



# Perform PCA ####
biotic_pca_result <- prcomp(biotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(biotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
biotic_pca_data <- as.data.frame(biotic_pca_result$x)
biotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
biotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(biotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
biotic_pe_biplot<- ggbiplot(biotic_pca_result, 
                            obs.scale = 1, 
                            var.scale = 2.5, 
                            groups = biotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                            ellipse = F, 
                            circle = TRUE,
                            varname.adjust = 1.5) +
  geom_mark_ellipse(aes(fill = biotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = biotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  scale_color_manual(values = c(c("PE477" = lighten("#0c1844", 0.4), "PE486" = lighten("#850000", 0.4))))+
  scale_fill_manual(values = c("PE477" = lighten("#0c1844", 0.6), "PE486" = lighten("#850000", 0.6))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-6, 6) +
  ylim(-8, 7)


biotic_pe_biplot

ggsave(biotic_pe_biplot, path = "./results/figures/", filename = "biotic_pe_vp_biplot.svg", dpi = 800, width = 9, height = 9)



#PCA

cols_to_remove <- c("Temperature", "Salinity", "Silicate", "Nitrate", "Nitrite", "Phosphate"
  )

a_biotic_numeric_columns_filtered <- numeric_columns_filtered %>%
  dplyr::select(-all_of(cols_to_remove), 
                -matches("Growth"), 
                -matches("decay"), 
                -matches("abs"),
                -matches("V1"), 
                -matches("V2"), 
                -matches("V3"))


# Standardize the data
a_biotic_numeric_columns_scaled <- scale(a_biotic_numeric_columns_filtered)



# Perform PCA ####
a_biotic_pca_result <- prcomp(a_biotic_numeric_columns_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA result
summary(a_biotic_pca_result)

# Scree plot to visualize the variance explained by each principal component
screeplot(pca_result, type = "lines")

# Create a data frame with PCA results
a_biotic_pca_data <- as.data.frame(a_biotic_pca_result$x)
a_biotic_pca_data$Station <- vp_pe_df$Station  # Add Station to the PCA data
a_biotic_pca_data$Cruise <- vp_pe_df$Cruise 

# Visualize the first two principal components
ggplot(a_biotic_pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Environmental Data", x = "Principal Component 1", y = "Principal Component 2")

# Create a biplot with arrows
a_biotic_pe_biplot<- ggbiplot(a_biotic_pca_result, 
                            obs.scale = 1, 
                            var.scale = 2.5, 
                            groups = a_biotic_pca_data$Cruise, # Assuming you want to color points by Cruise
                            ellipse = F, 
                            circle = F,
                            varname.adjust = 1.5) +
  # geom_mark_ellipse(aes(fill = a_biotic_pca_data$Cruise, color = NA), alpha = 0.1)  +
  geom_text_repel(aes(x = PC1, y = PC2, label = Station), data = a_biotic_pca_data, size = 5, vjust = -0.5,
                  max.overlaps = Inf)+
  scale_color_manual(values = c(c("PE477" ="#0c1844", "PE486" = "red")))+
  scale_fill_manual(values = c("PE477" = "#0c1844", "PE486" = "red")) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(#title = "PCA Biplot",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +
  xlim(-6, 6) +
  ylim(-8, 7)


a_biotic_pe_biplot
