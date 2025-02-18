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

pe_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7)


variables_to_test <- c("Temperature", "Salinity",
                       "Density", "Conductivity", "Turbidity",
                       "Nitrate", "Phosphate", "Silicate",
                       "Oxygen", "Chlorophyll",
                       "Total_Bacteria", 
                       "HNA", "LNA", "Cyanobacteria",
                       "Total_Viruses", "V1", "V2", "V3",
                       "VBR", "HNALNA")



# test

# Load necessary libraries
library(ggplot2)
library(ggpmisc)  # For regression equation

ggplot(pe_df_7m, aes(x = Temperature, y = Corrected_VP_Lytic)) +
geom_point(alpha = 0.6, color = "blue") +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Regression line
  stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), 
                                 after_stat(p.value), sep = "~~~")), 
               formula = y ~ x, 
               parse = TRUE, size = 5, color = "black", 
               label.x.npc = "right", label.y.npc = "top") +  # Position the text
  theme_minimal() +
  labs(title = "Corrected VP Lytic vs Temperature (Complete Data)",
       x = "Temperature (°C)", 
       y = "Corrected VP Lytic") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))


# 2.0 Lytic viral production ####
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)  # For pivoting data
library(ggpubr)
library(ggpmisc)  # For regression equation

# Filter datasets
pe_df_complete <- pe_df_7m  # Full dataset
pe_df_PE477 <- pe_df_7m %>% filter(Location == "PE477")
pe_df_PE486 <- pe_df_7m %>% filter(Location == "PE486")

# Function to create plots for different datasets
plot_faceted <- function(df, title) {
  # Convert to long format
  df_long <- df %>%
    select(Corrected_VP_Lytic, all_of(variables_to_test)) %>%
    pivot_longer(cols = -Corrected_VP_Lytic, names_to = "Variable", values_to = "Value")
  
  # Plot
  ggplot(df_long, aes(x = Value, y = Corrected_VP_Lytic)) +
    geom_point(alpha = 0.6, color = "blue") +  # Scatter points
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Regression line
    facet_wrap(~ Variable, scales = "free_x") +  # Faceting by variable
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   "\nR² = ", after_stat(rr.label),
                                   "\np = ", after_stat(p.value))), 
                 formula = y ~ x, 
                 parse = TRUE, size = 3, color = "black") +  # Equation, R², p-value
    theme_minimal() +
    labs(title = title, x = "Value", y = "Corrected VP Lytic") +
    theme(strip.text = element_text(size = 12, face = "bold"))
}

# Generate and display plots
plot_complete <- plot_faceted(pe_df_complete, "Corrected VP Lytic vs Viral Variables (Complete Data)")
plot_PE477 <- plot_faceted(pe_df_PE477, "Corrected VP Lytic vs Viral Variables (Location: PE477)")
plot_PE486 <- plot_faceted(pe_df_PE486, "Corrected VP Lytic vs Viral Variables (Location: PE486)")

# Display plots
print(plot_complete)
print(plot_PE477)
print(plot_PE486)



#########


# Load necessary libraries
library(ggplot2)
library(ggpmisc)
library(tidyr)
library(dplyr)

# Define datasets
pe_df_complete <- pe_df_7m  # Full dataset
pe_df_PE477 <- pe_df_7m %>% filter(Location == "PE477")  # PE477 only
pe_df_PE486 <- pe_df_7m %>% filter(Location == "PE486")  # PE486 only

# Add a new column to identify dataset
pe_df_complete$Dataset <- "Complete"
pe_df_PE477$Dataset <- "PE477"
pe_df_PE486$Dataset <- "PE486"

# Combine all datasets
pe_df_combined <- bind_rows(pe_df_complete, pe_df_PE477, pe_df_PE486)

# Convert the dataset to long format for faceting
pe_df_long <- pe_df_combined %>%
  dplyr::filter(Corrected_VP_Lytic > 0) %>%
  select(Dataset, Corrected_VP_Lytic, all_of(variables_to_test)) %>%
  pivot_longer(cols = -c(Corrected_VP_Lytic, Dataset), names_to = "Variable", values_to = "Value")

# Loop through each variable and save plots separately
for (var in variables_to_test) {
  plot_data <- pe_df_long %>% filter(Variable == var)  # Filter for the current variable
  
  p <- ggplot(plot_data, aes(x = Value, y = Corrected_VP_Lytic)) +
    geom_point(alpha = 0.6, color = "blue") +  # Scatter points
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Regression line
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(rr.label), 
                                   after_stat(p.value), sep = "~~~")), 
                 formula = y ~ x, 
                 parse = TRUE, size = 3) +  # Regression stats
    facet_grid(. ~ Dataset, scales = "free_x") +  # 3 columns, one per dataset
    theme_test() +
    labs(title = paste("Corrected VP Lytic vs", var),
         x = var, 
         y = "Corrected VP Lytic") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.background = element_rect(fill = "white"))
  
  # Save the plot
  ggsave(filename = paste0("./figures/LM/Lytic/VP_Lytic_vs_", var, ".png"), plot = p, width = 10, height = 5, dpi = 300)
  
  # Print confirmation
  print(paste("Saved:", paste0("VP_Lytic_vs_", var, ".png")))
}



# Convert the dataset to long format for faceting
pe_df_long <- pe_df_combined %>%
  dplyr::filter(Corrected_VP_Lysogenic>0) %>%
  select(Dataset, Corrected_VP_Lysogenic, all_of(variables_to_test)) %>%
  pivot_longer(cols = -c(Corrected_VP_Lysogenic, Dataset), names_to = "Variable", values_to = "Value")

# Loop through each variable and save plots separately
for (var in variables_to_test) {
  plot_data <- pe_df_long %>% filter(Variable == var)  # Filter for the current variable
  
  p <- ggplot(plot_data, aes(x = Value, y = Corrected_VP_Lysogenic)) +
    geom_point(alpha = 0.6, color = "blue") +  # Scatter points
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Regression line
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(rr.label), 
                                   after_stat(p.value), sep = "~~~")), 
                 formula = y ~ x, 
                 parse = TRUE, size = 3) +  # Regression stats
    facet_grid(. ~ Dataset, scales = "free_x") +  # 3 columns, one per dataset
    theme_test() +
    labs(title = paste("Corrected VP Lysogenic vs", var),
         x = var, 
         y = "Corrected VP Lysogenic") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.background = element_rect(fill = "white"))
  
  # Save the plot
  ggsave(filename = paste0("./figures/LM/Lysogenic/VP_Lysogenic_vs_", var, ".png"), plot = p, width = 10, height = 5, dpi = 300)
  
  # Print confirmation
  print(paste("Saved:", paste0("VP_Lysogenic_vs_", var, ".png")))
}



# Convert the dataset to long format for faceting
pe_df_long <- pe_df_combined %>%
  select(Dataset, decay_rate_linear , all_of(variables_to_test)) %>%
  pivot_longer(cols = -c(decay_rate_linear , Dataset), names_to = "Variable", values_to = "Value")

# Loop through each variable and save plots separately
for (var in variables_to_test) {
  plot_data <- pe_df_long %>% filter(Variable == var)  # Filter for the current variable
  
  p <- ggplot(plot_data, aes(x = Value, y = decay_rate_linear )) +
    geom_point(alpha = 0.6, color = "blue") +  # Scatter points
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Regression line
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(rr.label), 
                                   after_stat(p.value), sep = "~~~")), 
                 formula = y ~ x, 
                 parse = TRUE, size = 3) +  # Regression stats
    facet_grid(. ~ Dataset, scales = "free_x") +  # 3 columns, one per dataset
    theme_test() +
    labs(title = paste("decay_rate_linear  vs", var),
         x = var, 
         y = "decay_rate_linear") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.background = element_rect(fill = "white"))
  
  # Save the plot
  ggsave(filename = paste0("./figures/LM/Decay/decay_rate_vs_", var, ".png"), plot = p, width = 10, height = 5, dpi = 300)
  
  # Print confirmation
  print(paste("Saved:", paste0("decay_rate_vs_", var, ".png")))
}



# Wiht outliers removed ####
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(tidyr)
library(broom)
library(patchwork)

#Define full dataset and subsets
pe_df_complete <- pe_df_7m  # Full dataset
pe_df_PE477 <- pe_df_7m %>% filter(Location == "PE477")
pe_df_PE486 <- pe_df_7m %>% filter(Location == "PE486")

# Label datasets
pe_df_complete$Dataset <- "Complete"
pe_df_PE477$Dataset <- "PE477"
pe_df_PE486$Dataset <- "PE486"

# Combine all datasets
pe_df_combined <- bind_rows(pe_df_complete, pe_df_PE477, pe_df_PE486)
# 

# Define predictor and response variables
predictor_vars <- c("Temperature", "Salinity", "Turbidity",
                    "Nitrate", "Phosphate", "Silicate",
                    "Total_Bacteria", "Total_Viruses",
                    "HNA", "LNA", "V1", "V2", "V3", "VBR", "HNALNA",
                    "Cyanobacteria", "Chlorophyll")  # Example predictors
response_vars <- c(
  "Corrected_VP_Lytic", "Corrected_VP_Lysogenic", "decay_rate_linear",
                   "Total_Bacteria", "Total_Viruses",
                   "HNA", "LNA", "V1", "V2", "V3", "VBR", "HNALNA",
                   "Cyanobacteria", "Chlorophyll",
                   "percent_bacterial_loss_day_burst_50", "percent_lysogeny_day_burst_50", "percent_decay_day_linear")  # Example responses

# Define variables where zero values should be removed
exclude_zero_vars <- c("Corrected_VP_Lytic", "Corrected_VP_Lysogenic", "decay_rate_linear")

# Base directory for saving plots
base_dir <- "./figures/LM/"

# Ensure base directory exists
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# Function to compute Cook’s Distance and Standardized Residuals
compute_outlier_metrics <- function(df, x_var, y_var) {
  model <- lm(as.formula(paste(y_var, "~", x_var)), data = df)  
  cooks_d <- cooks.distance(model)
  std_resid <- rstandard(model)  
  
  threshold_cooks <- 4 / nrow(df)  # Cook's Distance Threshold
  threshold_resid <- 2  # Standardized Residual Threshold
  
  df$Cooks_Distance <- cooks_d
  df$Std_Residuals <- std_resid
  df$Outlier <- as.character(ifelse((cooks_d > threshold_cooks) & (abs(std_resid) > threshold_resid), "Outlier", "Normal"))
  
  return(df)
}

# Function to remove outliers based on both Cook’s Distance & Standardized Residuals
remove_outliers_combined <- function(df, x_var, y_var) {
  df_outlier_metrics <- compute_outlier_metrics(df, x_var, y_var)
  df_cleaned <- df_outlier_metrics %>% dplyr::filter(Outlier == "Normal")
  return(df_cleaned)
}

# Loop through each response variable and predictor
for (y_var in response_vars) {
  
  # Subset data **only for the current response variable** from pe_df_combined
  df_filtered <- pe_df_combined
  if (y_var %in% exclude_zero_vars) {
    df_filtered <- df_filtered %>%
      dplyr::filter(.data[[y_var]] != 0) %>%
      dplyr::select(Dataset, all_of(predictor_vars), all_of(y_var))  # Retain `Dataset`
  }
  
  # Ensure `Dataset` column exists after filtering
  if (!"Dataset" %in% colnames(df_filtered)) {
    stop(paste("Error: `Dataset` column missing after filtering for", y_var, "- Check dataset structure."))
  }
  
  # Check if the dataset has been completely filtered out
  if (nrow(df_filtered) == 0) {
    message(paste("Skipping", y_var, "as all values were removed during filtering."))
    next  # Skip to the next response variable
  }
  
  # Create directory for response variable inside base directory
  response_dir <- paste0(base_dir, y_var, "/")
  dir.create(response_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (x_var in predictor_vars) {
    
    # Compute outliers and clean data
    pe_df_outliers <- df_filtered %>%
      dplyr::group_by(Dataset) %>%
      dplyr::group_modify(~ compute_outlier_metrics(.x, x_var, y_var))
    
    pe_df_cleaned <- df_filtered %>%
      dplyr::group_by(Dataset) %>%
      dplyr::group_modify(~ remove_outliers_combined(.x, x_var, y_var))
    
    # Ensure `Outliers` column is consistent across datasets
    pe_df_outliers <- pe_df_outliers %>% mutate(Outliers = "Original")
    pe_df_cleaned <- pe_df_cleaned %>% mutate(Outliers = "No Outliers")
    
    # Generate regression plots
    reg_orig_complete <- plot_regression(pe_df_outliers %>% dplyr::filter(Dataset == "Complete"), x_var, y_var, "Complete (Original)")
    reg_orig_PE477 <- plot_regression(pe_df_outliers %>% dplyr::filter(Dataset == "PE477"), x_var, y_var, "PE477 (Original)")
    reg_orig_PE486 <- plot_regression(pe_df_outliers %>% dplyr::filter(Dataset == "PE486"), x_var, y_var, "PE486 (Original)")
    
    reg_clean_complete <- plot_regression(pe_df_cleaned %>% dplyr::filter(Dataset == "Complete"), x_var, y_var, "Complete (Outliers Removed)")
    reg_clean_PE477 <- plot_regression(pe_df_cleaned %>% dplyr::filter(Dataset == "PE477"), x_var, y_var, "PE477 (Outliers Removed)")
    reg_clean_PE486 <- plot_regression(pe_df_cleaned %>% dplyr::filter(Dataset == "PE486"), x_var, y_var, "PE486 (Outliers Removed)")
    
    # Generate diagnostic plots
    outlier_plot_complete <- plot_outlier_metrics(pe_df_outliers %>% dplyr::filter(Dataset == "Complete"), "Complete", x_var, y_var)
    outlier_plot_PE477 <- plot_outlier_metrics(pe_df_outliers %>% dplyr::filter(Dataset == "PE477"), "PE477", x_var, y_var)
    outlier_plot_PE486 <- plot_outlier_metrics(pe_df_outliers %>% dplyr::filter(Dataset == "PE486"), "PE486", x_var, y_var)
    
    # Arrange plots in a 3-column x 4-row grid
    final_plot <- (reg_orig_complete | reg_orig_PE477 | reg_orig_PE486) /
      (reg_clean_complete | reg_clean_PE477 | reg_clean_PE486) /
      (outlier_plot_complete | outlier_plot_PE477 | outlier_plot_PE486)
    
    # Save plot in respective folder
    plot_path <- paste0(response_dir, y_var, "_vs_", x_var, ".png")
    ggsave(plot_path, final_plot, width = 20, height = 10, dpi = 300)
    
    print(paste("Saved:", plot_path))
  }
}


# Explorin lytic viral production arte and LNA relation ####
lytic_df <- pe_df_7m %>%
  dplyr::filter(Corrected_VP_Lytic > 0)

lna_lytic_model <- lm(Corrected_VP_Lytic ~ LNA, data = lytic_df)
par(mfrow = c(2, 2)) 
plot(lna_lytic_model)

summary(lna_lytic_model)

lna_lytic_model <- lm(Corrected_VP_Lytic ~ LNA, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 3, 12.2)))
par(mfrow = c(2, 2)) 
plot(lna_lytic_model)

summary(lna_lytic_model)


ggplot(data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 3, 12.2)),
       aes(x = LNA,
           y = Corrected_VP_Lytic)) +
  geom_point()



# exploring lytic viral lproduction arte with HNA


hna_lytic_model <- lm(Corrected_VP_Lytic ~ HNA, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 2, 12.2)))
par(mfrow = c(2, 2)) 
plot(hna_lytic_model)

summary(hna_lytic_model)


ggplot(data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 2, 12.2)),
       aes(x = LNA,
           y = Corrected_VP_Lytic)) +
  geom_point()


# exploring lytic viral lproduction arte with total bacteria


bacteria_lytic_model <- lm(Corrected_VP_Lytic ~ Total_Bacteria, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 3, 12.2)))
par(mfrow = c(2, 2)) 
plot(bacteria_lytic_model)

summary(bacteria_lytic_model)


ggplot(data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 2, 12.2)),
       aes(x = LNA,
           y = Corrected_VP_Lytic)) +
  geom_point()


# exploring lytic viral lproduction arte with total VIRUSES


viruses_lytic_model <- lm(Corrected_VP_Lytic ~ Total_Viruses, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 6, 12.2)))
par(mfrow = c(2, 2)) 
plot(viruses_lytic_model)

summary(viruses_lytic_model)


ggplot(data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 6, 12.2)),
       aes(x = Total_Viruses,
           y = Corrected_VP_Lytic)) +
  geom_point()



viruses_lytic_model <- lm(Corrected_VP_Lytic ~ VBR, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 3, 6, 9,  12.2)))
par(mfrow = c(2, 2)) 
plot(viruses_lytic_model)

summary(viruses_lytic_model)


ggplot(data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 3, 6, 9, 12.2)),
       aes(x = VBR,
           y = Corrected_VP_Lytic)) +
  geom_point()


model_robust <- rlm(Corrected_VP_Lytic ~ VBR, data = lytic_df %>% dplyr::filter(!Station_Number %in% c(1, 12.2)))
par(mfrow = c(2, 2)) 
plot(model_robust)

summary(model_robust)


model_poly <- rlm(Corrected_VP_Lytic ~ VBR + I(VBR^2), data = lytic_df %>% dplyr::filter(!Station_Number %in% c()))
summary(model_poly)
par(mfrow = c(2, 2)) 
plot(model_poly)


weights <- 1 / (abs(resid(model_poly)) + 1)
model_weighted <- rlm(Corrected_VP_Lytic ~ VBR + I(VBR^2), data = lytic_df %>% dplyr::filter(!Station_Number %in% c()), weights = weights)
summary(model_weighted)

par(mfrow=c(2,2))
plot(model_weighted)






lysogenic_df <- pe_df_7m %>% dplyr::filter(Corrected_VP_Lysogenic > 0)

model_poly <- rlm(Corrected_VP_Lysogenic ~ VBR + I(VBR^2) , data = lysogenic_df %>% dplyr::filter(!Station_Number %in% c()))
summary(model_poly)
par(mfrow = c(2, 2)) 
plot(model_poly)


weights <- 1 / (abs(resid(model_poly)) + 1)
model_weighted <- rlm(Corrected_VP_Lysogenic ~ VBR + I(VBR^2), data = lysogenic_df %>% dplyr::filter(!Station_Number %in% c(6)), weights = weights)
summary(model_weighted)

par(mfrow=c(2,2))
plot(model_weighted)


#####weigthed polynomail regression
# Load necessary libraries
library(ggplot2)
library(MASS)
library(car)

# Fit the initial polynomial regression with weights
weights <- 1 / (abs(lm(Corrected_VP_Lytic ~ poly(VBR, 2), data = lytic_df)$residuals) + 0.1) # Avoid zero division
weighted_poly_model <- lm(Corrected_VP_Lytic ~ poly(VBR, 2), data = lytic_df, weights = weights)

# Summary of the model
summary(weighted_poly_model)
ggplot(lytic_df, aes(x = VBR, y = Corrected_VP_Lytic)) +
  geom_point(color = "black") +  # Data points
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "blue") +
  labs(title = "Weighted Polynomial Regression (Initial)",
       x = "VBR",
       y = "Corrected VP Lytic") +
  theme_minimal()

influencePlot(weighted_poly_model)


# Remove outliers 3 & 8
df_filtered <- lytic_df[-c(3, 8), ]  

# Recalculate weights without the removed points
weights_filtered <- 1 / (abs(lm(Corrected_VP_Lytic ~ poly(VBR, 2), data = df_filtered)$residuals) + 0.1)

# Fit new weighted polynomial regression
weighted_poly_filtered <- lm(Corrected_VP_Lytic ~ poly(VBR, 2), data = df_filtered, weights = weights_filtered)

# Summary of the new model
summary(weighted_poly_filtered)

AIC(weighted_poly_model, weighted_poly_filtered)

ggplot(df_filtered, aes(x = VBR, y = Corrected_VP_Lytic)) +
  geom_point(color = "black") +  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "blue") +
  labs(title = "Refitted Weighted Polynomial Regression (Outliers Removed)",
       x = "VBR",
       y = "Corrected VP Lytic") +
  theme_minimal()
shapiro.test(resid(weighted_poly_filtered))  # Normality check
ncvTest(weighted_poly_filtered)  # Homoscedasticity check
influencePlot(weighted_poly_filtered)  # New influence plot




### Turbidity 

# Select relevant variables
cor_matrix <- cor(pe_df[, c("Turbidity", "Salinity", "Chlorophyll", "Temperature")], use = "complete.obs")

# Display correlation matrix
print(cor_matrix)

# Visualize correlation matrix
library(ggcorrplot)
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

# Fit regression model
turbidity_model <- lm(Turbidity ~ Salinity + Chlorophyll + Temperature, data = pe_df)

# Summary of the model
summary(turbidity_model)

# Fit regression model
turbidity_model <- lm(Turbidity ~ Chlorophyll , data = pe_df %>% dplyr::filter(!Station_Number %in% c(8)))
par(mfrow=c(2,2))
plot(turbidity_model)
# Summary of the model
summary(turbidity_model)


# Fit regression model
chlorophyll_model <- lm(Chlorophyll ~ Temperature * Turbidity, data = lytic_df %>% dplyr::filter(!Station_Number %in% c()))
par(mfrow=c(2,2))
plot(chlorophyll_model)
# Summary of the model
summary(chlorophyll_model)

