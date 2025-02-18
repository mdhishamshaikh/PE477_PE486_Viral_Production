

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv")

physicochemical_params <- c("Temperature", "Salinity", "Density", "Conductivity", 
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
               "percent_lysogeny_burst_50")

sample_params <- c("Location", "Station_Number", "Depth", 
                   "sample_tag", "Latitude", "Longitude",
                   "Season")


# Dispersion tests to see if the variance in a season is higher thanin another
# Prepare PCA data (numeric only)
dispersion_pca_data <- pe_df %>%
  dplyr::select(all_of(physicochemical_params)) %>%
  select(where(is.numeric))

# Extract season labels
season_variable <- pe_df$Season  # Replace 'Season' with your actual column name

# Perform Beta Dispersion Test
betadisper_test <- betadisper(dist(dispersion_pca_data), season_variable)
permutest(betadisper_test) 




pe_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7)



# Dispersion tests to see if the variance in a season is higher thanin another
# Prepare PCA data (numeric only)
dispersion_pca_data <- pe_df %>%
  dplyr::select(all_of(biological_params)) %>%
  select(where(is.numeric))

# Extract season labels
season_variable <- pe_df$Season  # Replace 'Season' with your actual column name

# Perform Beta Dispersion Test
betadisper_test <- betadisper(dist(dispersion_pca_data), season_variable)
permutest(betadisper_test) 




dispersion_pca_pc_data <- pe_df_7m %>%
  dplyr::select(all_of(c(physicochemical_params))) %>%
  select(where(is.numeric))

# Extract season labels
season_variable <- pe_df_7m$Season  # Replace 'Season' with your actual column name

# Perform Beta Dispersion Test
betadisper_test_pc <- betadisper(dist(dispersion_pca_pc_data), season_variable)
permutest(betadisper_test_pc) 



dispersion_pca_pc_data <- pca_pc_df %>%
  dplyr::select(-all_of(sample_params)) %>%
  select(where(is.numeric))

# Extract season labels
season_variable <- pca_pc_labels$Season  # Replace 'Season' with your actual column name

# Perform Beta Dispersion Test
betadisper_test_pc <- betadisper(dist(dispersion_pca_pc_data), season_variable)
permutest(betadisper_test_pc)



library(vegan)

# Define response (viral production) and explanatory variables
response <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  select(percent_bacterial_loss_day_burst_30, percent_lysogeny_burst_30)

explanatory <- pe_df %>%
  dplyr::filter(Depth == 7) %>%
  select(Temperature, Turbidity, Salinity, Conductivity, Density, Total_Bacteria, Total_Viruses)

# Run RDA
rda_model <- rda(response ~ ., data = explanatory)
summary(rda_model)

# Plot RDA results
plot(rda_model, scaling = 2)
summary(rda_model)$biplot


# Linear regression
# Lytic viral production model
lm_lytic <- lm(percent_bacterial_loss_day_burst_30 ~  Total_Bacteria + Total_Viruses + VBR, data = pe_df %>% dplyr::filter(Depth == 7))
summary(lm_lytic)

# Lysogenic viral production model
lm_lysogenic <- lm(percent_lysogeny_burst_30 ~ Temperature + Season + Turbidity + Salinity + Conductivity + Density + Total_Bacteria + Total_Viruses, data = pe_df %>% dplyr::filter(Depth == 7))
summary(lm_lysogenic)

lm_lytic <- lm(Corrected_VP_Lytic ~ Temperature+Turbidity, data = pe_df %>% dplyr::filter(Depth == 7))
summary(lm_lytic)
lm_lysogenic <- lm(Corrected_VP_Lysogenic ~ Temperature+Turbidity, data = pe_df %>% dplyr::filter(Depth == 7))
summary(lm_lysogenic)

# Non linear

library(mgcv)
gam_lytic <- gam(percent_bacterial_loss_day_burst_30 ~ s(Total_Bacteria) + s(Total_Viruses) + s(Turbidity) +s(VBR), data = pe_df %>% filter(Depth == 7))
summary(gam_lytic)
plot(gam_lytic)



# pe rseason
lm_lytic_pe477 <- lm(Corrected_VP_Lytic ~  Temperature, data = pe_df %>% dplyr::filter(Depth == 7,Location == "PE477"))
summary(lm_lytic_pe477)
lm_lytic_pe486 <- lm(Corrected_VP_Lytic ~  Temperature, data = pe_df %>% dplyr::filter(Depth == 7,Location == "PE486"))
summary(lm_lytic_pe486)


lm_lysogenic<- lm(Corrected_VP_Lysogenic ~  Temperature, data = pe_df %>% dplyr::filter(Depth == 7, Corrected_VP_Lysogenic > 0))
summary(lm_lysogenic)


library(MASS)

# Start with a full model (all potential predictors)
full_model <- lm(Corrected_VP_Lytic ~ ., data = pe_df_7m %>% dplyr::select(Corrected_VP_Lytic, Temperature, Turbidity, Salinity, Conductivity, Density, Total_Bacteria, Total_Viruses))

# Stepwise selection (both forward & backward)
best_model_aic <- stepAIC(full_model, direction = "both", trace = FALSE)

# Print selected model summary
summary(best_model)


library(leaps)

# Perform exhaustive search
subset_model <- regsubsets(Corrected_VP_Lytic ~ ., data = pe_df_7m %>% dplyr::select(Corrected_VP_Lytic, Temperature, Turbidity, Salinity, Conductivity, Density, Total_Bacteria, Total_Viruses), nvmax = 5)

# Print summary of best models
summary(subset_model)

summary(subset_model)$adjr2  # Adjusted R² values
summary(subset_model)$bic    # Bayesian Information Criterion (BIC)
summary(subset_model)$cp     # Mallow’s Cp


# Best subset selection

library(leaps)

# Perform exhaustive search
subset_model <- regsubsets(Corrected_VP_Lytic ~ ., data = pe_df_7m %>% 
                             dplyr::select(Corrected_VP_Lytic, 
                                           Temperature, Turbidity, 
                                           Salinity,  
                                           Density, Max_Depth,
                                           Nitrate, Phosphate, Silicate,
                                           Oxygen, Fluorescence,
                                           Total_Bacteria, Total_Viruses), nvmax = 5)

# Print summary of best models
summary(subset_model)
summary(subset_model)$adjr2  # Adjusted R² values
summary(subset_model)$bic    # Bayesian Information Criterion (BIC)
summary(subset_model)$cp     # Mallows' Cp

best_model <- lm(Corrected_VP_Lytic ~ Temperature + Turbidity + Nitrate, data = pe_df_7m)
summary(best_model)


ggplot(data = pe_df_7m %>% dplyr::filter(Corrected_VP_Lytic> 0),
       aes(x = Temperature,
           y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~ Season, scales = "free_x")

ggplot(data = pe_df_7m %>% dplyr::filter(Corrected_VP_Lytic> 0),
       aes(x = Turbidity,
           y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Season, scales = "free_x")

ggplot(data = pe_df_7m %>% dplyr::filter(Corrected_VP_Lytic> 0),
       aes(x = Temperature,
           y = Turbidity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Season, scales = "free_x")

summary(lm(Corrected_VP_Lytic ~ Turbidity + Temperature + Season + Season*Temperature +Season*Turbidity, data = pe_df_7m))

ggplot(data = pe_df_7m,
       aes(x = Max_Depth,
           y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~ Season, scales = "free_x")

ggplot(data = pe_df_7m %>% dplyr::filter(Corrected_VP_Lysogenic> 0),
       aes(x = Temperature,
           y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~ Season, scales = "free_x")

ggplot(data = pe_df_7m %>% dplyr::filter(Corrected_VP_Lysogenic> 0),
       aes(x = Turbidity,
           y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Season, scales = "free_x")



# Filter data for each season
autumn_data <- pe_df_7m %>% filter(Location == "PE477")
spring_data <- pe_df_7m %>% filter(Location == "PE486")

# Fit linear models separately for each season
lm_autumn <- lm(Corrected_VP_Lytic ~ Max_Depth , data = autumn_data)
lm_spring <- lm(Corrected_VP_Lytic ~ Max_Depth , data = spring_data)

# Print model summaries
summary(lm_autumn)
summary(lm_spring)

# Compare adjusted R² values
adj_r2_autumn <- summary(lm_autumn)$adj.r.squared
adj_r2_spring <- summary(lm_spring)$adj.r.squared

print(paste("Adjusted R² for Autumn:", adj_r2_autumn))
print(paste("Adjusted R² for Spring:", adj_r2_spring))

anova(lm_autumn, lm_spring)




