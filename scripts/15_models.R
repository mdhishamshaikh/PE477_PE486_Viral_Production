# AIM: To explain viral production rates #####

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

# Some housekeeping first
# I will replace zeros in nitrate, phosphate, and silicate to half of theri detection limits
dl_nitrate <- 0.1529
dl_phosphate <- 0.0152
dl_silicate <- 0.0711

pe_df <- pe_df %>%
  mutate(
    Nitrate = ifelse(Nitrate == 0, dl_nitrate/2, Nitrate),
    Silicate = ifelse(Silicate == 0, dl_phosphate/2, Silicate),
    Phosphate = ifelse(Phosphate == 0, dl_silicate/2, Phosphate)
  )

# Subsetting for 7 m as this is the only depth where we performed viral production assays
pe_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7)
  


# 2.0 Lytic Viral Production ####

# Filtering for stations with >0 production rates as these represent successful experiments
lytic_pe_df <- pe_df_7m %>%
  dplyr::filter(Corrected_VP_Lytic > 0)

# Understandng correlations between lytic viral production rates and other variables

lytic_pe_corr_df <- lytic_pe_df %>%
  dplyr::select(Corrected_VP_Lytic,
                Temperature, Salinity, Turbidity,
                Nitrate, Phosphate, Silicate,
                Chlorophyll, Total_Bacteria, Total_Viruses, VBR, Cyanobacteria) 

lytic_pe_corr_matrix <- lytic_pe_corr_df %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

lytic_pe_melted_corr <- reshape2::melt(lytic_pe_corr_matrix)

# Plot the heatmap
lytic_pe_corr_plot<- ggplot(data = lytic_pe_melted_corr, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation plot: Lytic viral production vs environemntal variables",
       x = NULL,
       y = NULL)
lytic_pe_corr_plot

# Finding highly correlates variable pairs
lytic_pe_highly_correlated_pairs <- lytic_pe_melted_corr %>%
  dplyr::filter(abs(value) > 0.95 & Var1 != Var2) %>%
  arrange(desc(abs(value)))   

print(lytic_pe_highly_correlated_pairs)
# Temeprature and turidity are highly correlated

# Creating full linear model

lytic_pe_full_model <- lm(Corrected_VP_Lytic ~ Temperature + Salinity + 
                   Turbidity +
                   Nitrate + Phosphate + Silicate + 
                     Chlorophyll + 
                     Total_Bacteria + Total_Viruses + VBR + 
                   Cyanobacteria, 
                 data = lytic_pe_df)

summary(lytic_pe_full_model)

# Variance Inflation factor
vif(lytic_pe_full_model)

# VIF failed
# Figuring out why
colSums(is.na(lytic_pe_corr_df)) # no NAs, which is good!
apply(lytic_pe_corr_df, 2, var) # The variance is not great for some variables but is not zero as such so thisis oky.

# Back to correlation plot.
# Removing turbidity
lytic_pe_full_model <- lm(Corrected_VP_Lytic ~ Temperature + Salinity + 
                            #Turbidity +
                            Nitrate + Phosphate + Silicate + 
                            Chlorophyll + 
                            Total_Bacteria + Total_Viruses + VBR + 
                            Cyanobacteria, 
                          data = lytic_pe_df)

summary(lytic_pe_full_model)

vif(lytic_pe_full_model)

# As VIF will change depending upon which variables we remove, here's a function that keeps calculating it until we have all variables with VIF < 10 removed
remove_high_vif <- function(model, data, threshold = 10) {
  vif_values <- vif(model)
  
  while (max(vif_values) > threshold) {
    worst_var <- names(vif_values)[which.max(vif_values)]  # Identifying highest VIF variable
    predictors <- names(coef(model))[-1]  # Excluding Intercept
    predictors <- setdiff(predictors, worst_var)  # Removing the worst VIF variable
    
    formula_str <- as.formula(paste("Corrected_VP_Lytic ~", paste(predictors, collapse = " + ")))
    model <- lm(formula_str, data = data)  # Refitting model without high-VIF variable
    vif_values <- vif(model)
  }
  
  return(model)
}

# Run the function
lytic_pe_final_model <- remove_high_vif(lytic_pe_full_model, lytic_pe_corr_df)
vif(lytic_pe_final_model)
summary(lytic_pe_final_model)

reduced_model <- lm(Corrected_VP_Lytic ~ Temperature + Salinity + 
                      Total_Bacteria + Total_Viruses +
                      as.factor(Season), 
                    data = lytic_pe_df)
summary(reduced_model)


diagnostic_plots_extended <- function(model) {
  par(mfrow = c(3, 2))  # 3x2 layout for 6 plots
  
  # 1. Residuals vs Fitted
  plot(model$fitted.values, model$residuals,
       xlab = "Fitted Values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red", lwd = 2)
  
  # 2. Histogram of Residuals
  hist(model$residuals, main = "Histogram of Residuals",
       xlab = "Residuals", col = "lightblue", breaks = 10)
  
  # 3. Q-Q Plot
  qqnorm(model$residuals, main = "Q-Q Plot of Residuals")
  qqline(model$residuals, col = "red", lwd = 2)
  
  # 4. Cook's Distance
  plot(cooks.distance(model), type = "h",
       main = "Cook's Distance", ylab = "Cook's Distance")
  abline(h = 1, col = "red")
  
  # 5. Leverage Plot
  plot(hatvalues(model), main = "Leverage Plot", ylab = "Leverage")
  abline(h = 2 * mean(hatvalues(model)), col = "red")
  
  # 6. Standardized Residuals vs Fitted
  plot(model$fitted.values, rstandard(model),
       xlab = "Fitted Values", ylab = "Standardized Residuals",
       main = "Standardized Residuals vs Fitted")
  abline(h = c(-2, 0, 2), col = "red", lwd = 2, lty = 2)
  
  par(mfrow = c(1, 1))  # Reset layout
}

# Run the extended version
diagnostic_plots_extended(lytic_pe_final_model)
diagnostic_plots_extended(reduced_model)

cooks_distances <- cooks.distance(reduced_model)
n <- nrow(lytic_pe_df)  # Number of observations
threshold <- 4/n  # Common threshold

outliers <- which(cooks_distances > threshold)  # Identify outlier indices
print(outliers)

cooks_distances <- cooks.distance(spring_model)
n <- nrow(lytic_pe_df)  # Number of observations
threshold <- 4/n  # Common threshold

outliers <- which(cooks_distances > threshold)  # Identify outlier indices
print(outliers)



######## Trying for smaller size
library(brms)  
# Load necessary libraries
library(ggplot2)
library(boot)  # For bootstrapping
# For Bayesian regression
library(dplyr)

# Use spring_df dataset
df <- spring_df  

### 1️⃣ SIMPLE LINEAR REGRESSION
lm_model <- lm(decay_rate_linear ~ Turbidity, data = df)
summary_lm <- summary(lm_model)
print(summary_lm)  # Print regression summary

### 2️⃣ BOOTSTRAPPED CONFIDENCE INTERVALS
boot_fn <- function(data, indices) {
  model <- lm(decay_rate_linear ~ Turbidity, data = data[indices, ])
  return(coef(model)[2])  # Return slope coefficient
}

set.seed(42)
boot_results <- boot(df, boot_fn, R = 1000)  # Perform 1000 bootstrap resamples
boot_ci <- boot.ci(boot_results, type = "bca")  # Compute bias-corrected confidence intervals
print(boot_ci)

### 3️⃣ BAYESIAN REGRESSION (Using brms)
bayesian_model <- brm(decay_rate_linear ~ Turbidity, data = df, chains = 4, iter = 2000)
summary(bayesian_model)  # Print Bayesian regression summary

### 📊 4️⃣ VISUALIZATION: Scatter Plot with Regression Line
ggplot(df, aes(x = Turbidity, y = decay_rate_linear)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Corrected VP Lytic vs Temperature (Spring Data)",
       x = "Temperature",
       y = "Corrected VP Lytic") +
  theme_minimal()





#
bacteria_pe_full_model <- lm(Total_Bacteria ~ Temperature *  as.factor(Season)  + Salinity + 
                            #Turbidity +
                           # Nitrate + Phosphate + Silicate + 
                           # Chlorophyll + 
                             Total_Viruses + VBR + 
                            Cyanobacteria +
                             as.factor(Season), 
                          data = pe_df_7m)
summary(bacteria_pe_full_model)
vif(bacteria_pe_full_model)


#
bacteria_pe_temp_model <- lm(Total_Bacteria ~ Temperature, 
                             data = autumn_df)
summary(bacteria_pe_temp_model)
vif(bacteria_pe_temp_model)


spring_df <- subset(pe_df_7m, Season == "Spring (Apr 2021)")
autumn_df <- subset(pe_df_7m, Season == "Autumn (Sept 2020)")

