# To explain biological variables over 3 depths #

source("scripts/0_source.R")

# 1.0 Importing combined data frame ####

pe_df <- read.csv("./results/PE477_PE486_3depths_combined.csv") %>%
  dplyr::mutate(Total_phyto = Total_phyto/1e+3,
                Total_Bacteria = Total_Bacteria/1e+6,
                Total_Viruses = Total_Viruses/1e+6,
                Cyanobacteria = Cyanobacteria/1e+3,
                HNA = HNA/1e+6,
                LNA = LNA/1e+6,
                V1 = V1/1e+6,
                V2 = V2/1e+6,
                V3 = V3/1e+6
  )


# Assigning vectors for parameters

physicochemical_params <- c("Temperature", "Salinity", #"Density", "Conductivity", 
                            "Turbidity", "Nitrate", "Phosphate", "Silicate",
                            "Max_Depth")

biological_params <- c(#"Oxygen", 
  "Chlorophyll", "Total_phyto",
  "Total_Bacteria", "HNA", "LNA", "Cyanobacteria",
  "Total_Viruses", "V1", "V2", "V3",
  "VBR", "HNALNA")


sample_params <- c("Location", "Station_Number", "Depth", 
                   "sample_tag", "Latitude", "Longitude",
                   "Season")

# Subsetting for only essential variables

pe_3d <- pe_df %>%
  dplyr::select(all_of(c(sample_params, physicochemical_params, biological_params))) %>%
  na.omit()

pe_3d_lr <- pe_3d %>%
  dplyr::select(-c("Location", "Station_Number", # "Depth", 
                   "sample_tag", "Latitude", "Longitude")) %>%
  dplyr::mutate(Season = as.factor(Season),
                Depth = as.factor(Depth))

is.na(pe_3d_lr)



# 2.0 Correlation ####


# Checking for correlation 
cor_matrix <- cor(pe_3d_lr[, sapply(pe_3d_lr, is.numeric)], 
                     use = "pairwise.complete.obs", method = "pearson")

# Melt the correlation matrix for heatmap visualization
melted_cor_matrix <- reshape2::melt(cor_matrix)

# Plot the correlation heatmap
cor_plot <- ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Heatmap: Physicochemical Variables",
       x = NULL, y = NULL)

# Display the heatmap
print(cor_plot)

highly_correlated_pairs <- melted_cor_matrix %>%
  dplyr::filter(abs(value) > 0.80 & Var1 != Var2) %>%
  arrange(desc(abs(value)))

# Print the pairs
print(highly_correlated_pairs)

# Removing phosphate, V1, V2, HNA, LNA, Cyanobacteria and temperature

low_corr_df <- pe_3d_lr %>%
  dplyr::select(-c(Season, Temperature, Phosphate, 
                   Cyanobacteria, HNA, LNA, V1, V2, V3, VBR, HNALNA))
# Checking for correlation 
low_cor_matrix <- cor(low_corr_df[, sapply(low_corr_df, is.numeric)], 
                  use = "pairwise.complete.obs", method = "pearson")

# Melt the correlation matrix for heatmap visualization
low_melted_cor_matrix <- reshape2::melt(low_cor_matrix)

ggplot(data = low_melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#003049", high = "#DC2828", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  labs(title = "Correlation Heatmap: Physicochemical Variables",
       x = NULL, y = NULL)
low_highly_correlated_pairs <- low_melted_cor_matrix %>%
  dplyr::filter(abs(value) > 0.80 & Var1 != Var2) %>%
  arrange(desc(abs(value)))

# Print the pairs
print(low_highly_correlated_pairs)




# 2.0 Stepwise linear regressionto explain total phytoplankton ####

phyto_full_model <- lm(Total_phyto ~ ., data=low_corr_df)
summary(phyto_full_model)

phyto_vif <- vif(phyto_full_model)
phyto_vif

# Forward selecction
phyto_null_model <- lm(Total_phyto ~ 1, data=low_corr_df)  # Start with an intercept-only model
phyto_full_model <- lm(Total_phyto ~ ., data=low_corr_df)  # Full model with all predictors

phyto_forward_model <- step(phyto_null_model, 
                               scope = list(lower = phyto_null_model, upper = phyto_full_model), 
                               direction = "forward")
summary(phyto_forward_model)

phyto_backward_model <- step(phyto_full_model, direction = "backward")
summary(phyto_backward_model)
phyto_stepwise_model <- step(phyto_full_model, direction = "both")
summary(phyto_stepwise_model)


par(mfrow=c(2,2))
plot(phyto_stepwise_model)


phyto_interaction_model <- lm(Total_phyto ~ Salinity + Turbidity * Total_Bacteria, data=low_corr_df)
summary(phyto_interaction_model)

par(mfrow=c(2,2))
plot(phyto_interaction_model)

phyto_log_model <- lm(log(Total_phyto) ~ Salinity + Turbidity * Total_Bacteria, data = low_corr_df)
summary(phyto_log_model)

par(mfrow=c(2,2))
plot(phyto_log_model)

AIC(phyto_stepwise_model, phyto_interaction_model, phyto_log_model)

anova(phyto_stepwise_model, phyto_interaction_model, phyto_log_model)
# Log excluded, iteraction term improves


# pERFORMING BOX COX TO VERIFY A LOG TRANSFORMATION WA SNEEDED
par(mfrow=c(1,1))
boxcox(lm(Total_phyto ~ Salinity + Turbidity * Total_Bacteria, data = low_corr_df), lambda = seq(-2, 2, by = 0.1))
# Lambda near zero.



# 3.0 Stepwise linear regression to explain total bacteria ####

bacteria_full_model <- lm(Total_Bacteria ~ ., data=low_corr_df)
summary(bacteria_full_model)

bacteria_vif <- vif(bacteria_full_model)
bacteria_vif

# Forward selecction
bacteria_null_model <- lm(Total_Bacteria ~ 1, data=low_corr_df)  # Start with an intercept-only model
bacteria_full_model <- lm(Total_Bacteria ~ ., data=low_corr_df)  # Full model with all predictors

bacteria_forward_model <- step(bacteria_null_model, 
                      scope = list(lower = bacteria_null_model, upper = bacteria_full_model), 
                      direction = "forward")
summary(bacteria_forward_model)

bacteria_backward_model <- step(bacteria_full_model, direction = "backward")
summary(bacteria_backward_model)
bacteria_stepwise_model <- step(bacteria_full_model, direction = "both")
summary(bacteria_stepwise_model)

bacteria_stepwise_model <- lm(Total_Bacteria ~ Turbidity + Nitrate + Total_phyto, data = low_corr_df)
summary(bacteria_stepwise_model)
par(mfrow=c(2,2))
plot(bacteria_stepwise_model)
vif(bacteria_stepwise_model)


bacteria_interaction_model <- lm(Total_Bacteria ~  Nitrate +  Turbidity * Total_phyto, data = low_corr_df)
summary(bacteria_interaction_model)
par(mfrow=c(2,2))
plot(bacteria_interaction_model)
vif(bacteria_interaction_model)

bacteria_poly_interaction_model <- lm(Total_Bacteria ~ Nitrate + Turbidity * poly(Total_phyto, 2), data = low_corr_df)
summary(bacteria_poly_interaction_model)
par(mfrow=c(2,2))
plot(bacteria_poly_interaction_model)
vif(bacteria_poly_interaction_model)
par(mfrow=c(1,1))
influencePlot(bacteria_poly_interaction_model)

AIC(bacteria_stepwise_model, bacteria_interaction_model, bacteria_poly_interaction_model)

anova(bacteria_stepwise_model, bacteria_interaction_model, bacteria_poly_interaction_model)
# Imporvement seen



# 4.0 Stepwise linear regression to explain total viruses ####

viruses_full_model <- lm(Total_Viruses ~ ., data=low_corr_df)
summary(viruses_full_model)

viruses_vif <- vif(viruses_full_model)
viruses_vif

# Forward selecction
viruses_null_model <- lm(Total_Viruses ~ 1, data=low_corr_df)  # Start with an intercept-only model
viruses_full_model <- lm(Total_Viruses ~ ., data=low_corr_df)  # Full model with all predictors

viruses_forward_model <- step(viruses_null_model, 
                               scope = list(lower = viruses_null_model, upper = viruses_full_model), 
                               direction = "forward")
summary(viruses_forward_model)

viruses_backward_model <- step(viruses_full_model, direction = "backward")
summary(viruses_backward_model)
viruses_stepwise_model <- step(viruses_full_model, direction = "both")
summary(viruses_stepwise_model)

viruses_stepwise_model <- lm(Total_Viruses ~  Max_Depth + Nitrate + Total_phyto, data = low_corr_df)
summary(viruses_stepwise_model)
par(mfrow=c(2,2))
plot(viruses_stepwise_model)
vif(viruses_stepwise_model)
bptest(viruses_stepwise_model)

viruses_log_model <- lm(log(Total_Viruses) ~ Max_Depth + Nitrate * Total_phyto, data = low_corr_df)
summary(viruses_log_model)

par(mfrow=c(2,2))
plot(viruses_log_model)
vif(viruses_log_model)
bptest(viruses_log_model)

AIC(viruses_stepwise_model, viruses_log_model )

# pERFORMING BOX COX TO VERIFY A LOG TRANSFORMATION WA SNEEDED
par(mfrow=c(1,1))
boxcox(lm(log(Total_Viruses) ~ Max_Depth + Nitrate * Total_phyto, data = low_corr_df), lambda = seq(-2, 2, by = 0.1))
# Lambda not near zero, log not needed





# 5.0 Model summaries ####

# All models
phyto_stepwise_model
phyto_interaction_model
phyto_log_model
bacteria_stepwise_model
bacteria_interaction_model
bacteria_poly_interaction_model
viruses_stepwise_model
viruses_log_model

# selected models 

# Total phytoplankton: 
phyto_log_model <- lm(log(Total_phyto) ~ Salinity + Turbidity * Total_Bacteria, data = low_corr_df)

# Total bacteria:
bacteria_poly_interaction_model <- lm(Total_Bacteria ~ Nitrate + Turbidity * poly(Total_phyto, 2), data = low_corr_df)

# Total viruses:
viruses_stepwise_model <- lm(Total_Viruses ~  Max_Depth + Nitrate + Total_phyto, data = low_corr_df)

# All combined

library(sjPlot)

tab_model(
  phyto_stepwise_model, phyto_interaction_model, phyto_log_model, 
  bacteria_stepwise_model, bacteria_interaction_model, bacteria_poly_interaction_model,
  viruses_stepwise_model, viruses_log_model,
  show.aic = TRUE, show.r2 = TRUE, 
  dv.labels = c("Phyto Stepwise", "Phyto Interaction", "Phyto Log",
                "Bacteria Stepwise", "Bacteria Interaction", "Bacteria Polynomial Interaction",
                "Viruses Stepwise", "Viruses Log"),
  title = "Model Summary: Phytoplankton, Bacteria, and Viruses"
)


# Phyto

tab_model(
  phyto_stepwise_model, phyto_interaction_model, phyto_log_model, 
  show.aic = TRUE, show.r2 = TRUE, 
  dv.labels = c("Phyto Stepwise", "Phyto Interaction", "Phyto Log"),
  title = "Model Summary: Phytoplankton"
)

# Bacteria

tab_model(
  bacteria_stepwise_model, bacteria_interaction_model, bacteria_poly_interaction_model,
  show.aic = TRUE, show.r2 = TRUE, 
  dv.labels = c("Bacteria Stepwise", "Bacteria Interaction", "Bacteria Polynomial Interaction"),
  title = "Model Summary: Bacteria"
)


# Viruses

tab_model(
  viruses_stepwise_model, viruses_log_model,
  show.ci = TRUE, # Add confidence intervals
  show.aic = TRUE, show.r2 = TRUE, 
  dv.labels = c("Viruses Stepwise", "Viruses Log"),
  title = "Model Summary: Viruses"
)























# 5.0 Stepwise linear regression to explain turbidity ####

turb_df <- cbind(low_corr_df, pe_3d_lr %>% na.omit() %>% dplyr::select(Season)) %>%
  dplyr::mutate(Season = as.factor(Season))


turbidity_full_model <- lm(Turbidity ~ ., data=turb_df)
summary(turbidity_full_model)

turbidity_vif <- vif(turbidity_full_model)
turbidity_vif

# Forward selecction
turbidity_null_model <- lm(Turbidity ~ 1, data=turb_df)  # Start with an intercept-only model
turbidity_full_model <- lm(Turbidity ~ ., data=turb_df)  # Full model with all predictors

turbidity_forward_model <- step(turbidity_null_model, 
                              scope = list(lower = turbidity_null_model, upper = turbidity_full_model), 
                              direction = "forward")
summary(turbidity_forward_model)

turbidity_backward_model <- step(turbidity_full_model, direction = "backward")
summary(turbidity_backward_model)
turbidity_stepwise_model <- step(turbidity_full_model, direction = "both")
summary(turbidity_stepwise_model)

turbidity_stepwise_model <- lm(Turbidity ~   Max_Depth + Nitrate + Season, data = turb_df)
summary(turbidity_stepwise_model)
par(mfrow=c(2,2))
plot(turbidity_stepwise_model)
vif(turbidity_stepwise_model)
bptest(turbidity_stepwise_model)

turbidity_interaction_model <- lm(log(Turbidity) ~ (Max_Depth + Nitrate) * Season, data = turb_df)
summary(turbidity_interaction_model)

par(mfrow=c(2,2))
plot(turbidity_interaction_model)
vif(turbidity_interaction_model)
bptest(turbidity_interaction_model)

AIC(turbidity_stepwise_model, turbidity_interaction_model )

