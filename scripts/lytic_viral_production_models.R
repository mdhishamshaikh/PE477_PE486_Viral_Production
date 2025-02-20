# To explain lytic viral production rates #

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

lytic_df <- pe_df_7m %>%
  dplyr::filter(Corrected_VP_Lytic > 0)

# 1.0 Lytic viral production rate vs Total bacteria ####

bacteria_lytic_lm <- lm(Corrected_VP_Lytic ~ Total_Bacteria, data = lytic_df)

par(mfrow=c(2,2))
plot(bacteria_lytic_lm)

summary(bacteria_lytic_lm)

weights <- 1 / (abs(bacteria_lytic_lm$residuals) + 0.1)  # Avoid division by zero

bacteria_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ Total_Bacteria, 
                                 data = lytic_df, 
                                 weights = weights)

par(mfrow=c(2,2))
plot(bacteria_lytic_weighted_lm)
summary(bacteria_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(bacteria_lytic_weighted_lm)

AIC(bacteria_lytic_lm, bacteria_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(bacteria_lytic_lm)
summary_weighted <- summary(bacteria_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(bacteria_lytic_lm)[1], 2), 
                        " + ", round(coef(bacteria_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(bacteria_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(bacteria_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = Total_Bacteria, y = Corrected_VP_Lytic, shape = Season)) +
  geom_point(size = 4) +
  scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$Total_Bacteria), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$Total_Bacteria), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
       x = "Total Bacteria", y = "Corrected VP Lytic") +
  theme_test()


# 2.0 Lytic viral production rate vs HNA ####

hna_lytic_lm <- lm(Corrected_VP_Lytic ~ HNA, data = lytic_df)

par(mfrow=c(2,2))
plot(hna_lytic_lm)

summary(hna_lytic_lm)

weights <- 1 / (abs(hna_lytic_lm$residuals) + 0.1)  # Avoid division by zero

hna_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ HNA, 
                                 data = lytic_df, 
                                 weights = weights)

par(mfrow=c(2,2))
plot(hna_lytic_weighted_lm)
summary(hna_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(hna_lytic_weighted_lm)

AIC(hna_lytic_lm, hna_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(hna_lytic_lm)
summary_weighted <- summary(hna_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(hna_lytic_lm)[1], 2), 
                        " + ", round(coef(hna_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(hna_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(hna_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = HNA, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$HNA), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$HNA), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "HNA bacteria", y = "Corrected VP Lytic") +
  theme_test()

# 3.0 Lytic viral production rate vs LNA ####

lna_lytic_lm <- lm(Corrected_VP_Lytic ~ LNA, data = lytic_df)

par(mfrow=c(2,2))
plot(lna_lytic_lm)

summary(lna_lytic_lm)

weights <- 1 / (abs(lna_lytic_lm$residuals) + 0.1)  # Avoid division by zero

lna_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ LNA, 
                            data = lytic_df, 
                            weights = weights)

par(mfrow=c(2,2))
plot(lna_lytic_weighted_lm)
summary(lna_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(lna_lytic_weighted_lm)

AIC(lna_lytic_lm, lna_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(lna_lytic_lm)
summary_weighted <- summary(lna_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(lna_lytic_lm)[1], 2), 
                        " + ", round(coef(lna_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(lna_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(lna_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = LNA, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$LNA), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$LNA), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "LNA bacteria", y = "Corrected VP Lytic") +
  theme_test()

# 4.0 Lytic viral production rate vs Total viruses ####

viruses_lytic_lm <- lm(Corrected_VP_Lytic ~ Total_Viruses, data = lytic_df)

par(mfrow=c(2,2))
plot(viruses_lytic_lm)

summary(viruses_lytic_lm)

weights <- 1 / (abs(viruses_lytic_lm$residuals) + 0.1)  # Avoid division by zero

viruses_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ Total_Viruses, 
                                 data = lytic_df, 
                                 weights = weights)

par(mfrow=c(2,2))
plot(viruses_lytic_weighted_lm)
summary(viruses_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(viruses_lytic_weighted_lm)

AIC(viruses_lytic_lm, viruses_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(viruses_lytic_lm)
summary_weighted <- summary(viruses_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(viruses_lytic_lm)[1], 2), 
                        " + ", round(coef(viruses_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(viruses_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(viruses_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = Total_Viruses, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$Total_Viruses), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$Total_Viruses), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Total Viruses", y = "Corrected VP Lytic") +
  theme_test()

# 5.0 Lytic viral production rate vs V1 ####

v1_lytic_lm <- lm(Corrected_VP_Lytic ~ V1, data = lytic_df)

par(mfrow=c(2,2))
plot(v1_lytic_lm)

summary(v1_lytic_lm)

weights <- 1 / (abs(v1_lytic_lm$residuals) + 0.1)  # Avoid division by zero

v1_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ V1, 
                                data = lytic_df, 
                                weights = weights)

par(mfrow=c(2,2))
plot(v1_lytic_weighted_lm)
summary(v1_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v1_lytic_weighted_lm)

AIC(v1_lytic_lm, v1_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(v1_lytic_lm)
summary_weighted <- summary(v1_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v1_lytic_lm)[1], 2), 
                        " + ", round(coef(v1_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v1_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(v1_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = V1, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$V1), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$V1), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V1", y = "Corrected VP Lytic") +
  theme_test()

# 6.0 Lytic viral production rate vs V2 ####

v2_lytic_lm <- lm(Corrected_VP_Lytic ~ V2, data = lytic_df)

par(mfrow=c(2,2))
plot(v2_lytic_lm)

summary(v2_lytic_lm)

weights <- 1 / (abs(v2_lytic_lm$residuals) + 0.1)  # Avoid division by zero

v2_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ V2, 
                           data = lytic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(v2_lytic_weighted_lm)
summary(v2_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v2_lytic_weighted_lm)

AIC(v2_lytic_lm, v2_lytic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(v2_lytic_lm)
summary_weighted <- summary(v2_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v2_lytic_lm)[1], 2), 
                        " + ", round(coef(v2_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v2_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(v2_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = V2, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$V2), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$V2), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V2", y = "Corrected VP Lytic") +
  theme_test()

# 7.0 Lytic viral production rate vs V3 ####

v3_lytic_lm <- lm(Corrected_VP_Lytic ~ V3, data = lytic_df)

par(mfrow=c(2,2))
plot(v3_lytic_lm)

summary(v3_lytic_lm)

weights <- 1 / (abs(v3_lytic_lm$residuals) + 0.1)  # Avoid division by zero

v3_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ V3, 
                           data = lytic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(v3_lytic_weighted_lm)
summary(v3_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v3_lytic_weighted_lm)

AIC(v3_lytic_lm, v3_lytic_weighted_lm)



# Get model summaries
summary_unweighted <- summary(v3_lytic_lm)
summary_weighted <- summary(v3_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v3_lytic_lm)[1], 2), 
                        " + ", round(coef(v3_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v3_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(v3_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = V3, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$V3), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$V3), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V3", y = "Corrected VP Lytic") +
  theme_test()

# 8.0 Lytic viral production rate vs Cyanobacteria ####
cyano_df_filtered <- lytic_df %>% dplyr::filter(!Station_Number %in% c(5, 6))

cyano_lytic_lm <- lm(Corrected_VP_Lytic ~ Cyanobacteria, data = cyano_df_filtered)

par(mfrow=c(2,2))
plot(cyano_lytic_lm)

summary(cyano_lytic_lm)

weights <- 1 / (abs(cyano_lytic_lm$residuals) + 0.1)  # Avoid division by zero

cyano_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ Cyanobacteria, 
                            data = cyano_df_filtered, 
                            weights = weights)

par(mfrow=c(2,2))
plot(cyano_lytic_weighted_lm)
summary(cyano_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(cyano_lytic_weighted_lm)

AIC(cyano_lytic_lm, cyano_lytic_weighted_lm)


# Get model summaries
summary_unweighted <- summary(cyano_lytic_lm)
summary_weighted <- summary(cyano_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(cyano_lytic_lm)[1], 2), 
                        " + ", round(coef(cyano_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(cyano_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(cyano_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(cyano_df_filtered, aes(x = Cyanobacteria, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(cyano_df_filtered$Cyanobacteria), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(cyano_df_filtered$Cyanobacteria), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Cyanobacteria", y = "Corrected VP Lytic") +
  theme_test()

# 9.0 Lytic viral production with VBR ####



vbr_lytic_lm <- lm(Corrected_VP_Lytic ~ VBR, data = lytic_df)

par(mfrow=c(2,2))
plot(vbr_lytic_lm)

summary(vbr_lytic_lm)

weights <- 1 / (abs(vbr_lytic_lm$residuals) + 0.1)  # Avoid division by zero

vbr_lytic_weighted_lm <- lm(Corrected_VP_Lytic ~ VBR, 
                           data = lytic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(vbr_lytic_weighted_lm)
summary(vbr_lytic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(vbr_lytic_weighted_lm)

AIC(vbr_lytic_lm, vbr_lytic_weighted_lm)
# weighted regression is better.  

# Get model summaries
summary_unweighted <- summary(vbr_lytic_lm)
summary_weighted <- summary(vbr_lytic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(vbr_lytic_lm)[1], 2), 
                        " + ", round(coef(vbr_lytic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(vbr_lytic_weighted_lm)[1], 2), 
                      " + ", round(coef(vbr_lytic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lytic_df, aes(x = VBR, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$VBR), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$VBR), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "VBR", y = "Corrected VP Lytic") +
  theme_test()


#POLYNOMIAL#
# Load required libraries
library(MASS)  # For robust regression
library(ggplot2)
library(car)   # For influence plots

# Fit an unweighted second-degree polynomial regression (LM)
vbr_lytic_poly <- lm(Corrected_VP_Lytic ~ poly(VBR, 2), data = lytic_df)

# Compute weights: Inverse of absolute residuals + small constant to avoid division by zero
weights <- 1 / (abs(vbr_lytic_poly$residuals) + 0.1)

# Fit weighted polynomial regression (LM)
vbr_lytic_weighted_poly <- lm(Corrected_VP_Lytic ~ poly(VBR, 2), 
                              data = lytic_df, 
                              weights = weights)

# Fit robust polynomial regression (RLM)
vbr_lytic_robust_poly <- rlm(Corrected_VP_Lytic ~ poly(VBR, 2), data = lytic_df)

# Compute new weights for robust regression
weights_robust <- 1 / (abs(resid(vbr_lytic_robust_poly)) + 0.1)

# Fit weighted robust polynomial regression (RLM with weights)
vbr_lytic_weighted_robust_poly <- rlm(Corrected_VP_Lytic ~ poly(VBR, 2), 
                                      data = lytic_df, 
                                      weights = weights_robust)

# Summaries
summary_unweighted <- summary(vbr_lytic_poly)
summary_weighted <- summary(vbr_lytic_weighted_poly)
summary_robust <- summary(vbr_lytic_robust_poly)
summary_weighted_robust <- summary(vbr_lytic_weighted_robust_poly)
summary_unweighted
summary_weighted
summary_robust
summary_weighted_robust
# Compare AIC values
AIC(vbr_lytic_poly, vbr_lytic_weighted_poly, vbr_lytic_robust_poly, vbr_lytic_weighted_robust_poly)

# Diagnostic Plots
par(mfrow=c(2,2))
plot(vbr_lytic_poly)
plot(vbr_lytic_weighted_poly)
plot(vbr_lytic_robust_poly)
plot(vbr_lytic_weighted_robust_poly)

# Influence plots
par(mfrow=c(1,1))
influencePlot(vbr_lytic_weighted_poly)
influencePlot(vbr_lytic_weighted_robust_poly)

# Extract model equations and statistics
eq_unweighted <- paste0("y = ", round(coef(vbr_lytic_poly)[1], 2), 
                        " + ", round(coef(vbr_lytic_poly)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(vbr_lytic_weighted_poly)[1], 2), 
                      " + ", round(coef(vbr_lytic_weighted_poly)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

eq_robust <- paste0("y = ", round(coef(vbr_lytic_robust_poly)[1], 2), 
                    " + ", round(coef(vbr_lytic_robust_poly)[2], 2), "*x\n",
                    "Adj R² = NA (RLM)", 
                    ", p = NA (RLM)")

eq_weighted_robust <- paste0("y = ", round(coef(vbr_lytic_weighted_robust_poly)[1], 2), 
                             " + ", round(coef(vbr_lytic_weighted_robust_poly)[2], 2), "*x\n",
                             "Adj R² = NA (RLM)", 
                             ", p = NA (RLM)")

# Create scatterplots with polynomial regression lines
ggplot(lytic_df, aes(x = VBR, y = Corrected_VP_Lytic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lytic_df$VBR), y = 8, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lytic_df$VBR), y = 7, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = "#008500", se = FALSE) +
  annotate("text", x =35, y = 8, 
           label = eq_robust, hjust = 0, size = 5, color = "#008500") +
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = "#856300", se = FALSE, 
              aes(weight = weights_robust)) +
  annotate("text", x = 35, y = 7, 
           label = eq_weighted_robust, hjust = 0, size = 5, color = "#856300") +
  labs(#title = "Polynomial Regression: Unweighted vs. Weighted vs. Robust",
    x = "VBR", y = "Corrected VP Lytic") +
  theme_test()


confint(vbr_lytic_poly)
confint(vbr_lytic_weighted_poly)
confint(vbr_lytic_robust_poly)  
confint(vbr_lytic_weighted_robust_poly)  

# 10.0 Facet plots for lytic viral pproduction WLS ####


long_lytic_df <- lytic_df %>%
  pivot_longer(cols = c(Total_Bacteria, HNA, LNA, Total_Viruses, V1, V2),
               names_to = "Predictor",
               values_to = "Predictor_Value") %>%
  dplyr::mutate(Predictor = factor(Predictor, levels = c("Total_Bacteria","HNA", "LNA",
                                                         "Total_Viruses", "V1", "V2"))) 
long_lytic_df <- long_lytic_df %>%
  group_by(Predictor) %>%
  mutate(
    # Fit an initial OLS model per predictor
    ols_model = list(lm(Corrected_VP_Lytic ~ Predictor_Value, data = cur_data())),
    
    # Extract residuals from the OLS model
    residuals_ols = abs(residuals(ols_model[[1]])) + 0.1,  # Adding small constant to avoid zero division
    
    # Compute weights as inverse of residuals
    weights = 1 / residuals_ols
  ) %>%
  ungroup() %>%
  dplyr::select(-ols_model, -residuals_ols)  # Remove intermediate columns

# Check if weights were added successfully
glimpse(long_lytic_df)

# Generate faceted WLS plot
lytic_wls_plot<- ggplot(long_lytic_df, aes(x = Predictor_Value, y = Corrected_VP_Lytic, shape = Season, fill = weights)) +
  # WLS Regression (Single Line per Facet)
  geom_smooth(method = "lm", formula = y ~ x, color = NA, se = T, size = 1, 
              aes(weight = weights, group = 1)) +
  
  geom_point(size = 3, stroke = 0.8, color = "black") +  # Ensures points have a black border
  scale_shape_manual(values = c(21, 24)) +  # Circle (21) & Triangle (24) for season
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Weight") +  # Shared weight color scale
  
  # WLS Regression (Single Line per Facet)
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = F, size = 1, 
              aes(weight = weights, group = 1)) +
  
  # Facet by predictor variable
  facet_wrap(~ Predictor, scales = "free_x", strip.position = "bottom",
             labeller = as_labeller(variable_labels, label_parsed)) +  
  
  # Labels (no title)
  labs(x = NULL, 
       y = expression("Lytic viral production rate " ~ (10^{5} ~ viruses ~ mL^{-1} ~ h^{-1})), 
       fill = "Weight") +
  
  # Theme adjustments for clean, publication-ready visuals
  theme_test(base_size = 15) +
  theme(legend.position = "right",
       panel.spacing.y = unit(0.5, "cm"),
      #panel.spacing.x = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 1.5),
        axis.line = element_line(size = 0), 
      axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.placement = "outside", 
      strip.switch.pad.wrap = unit(-0.1, "cm"),
      strip.text = element_text(color = "black"),
      aspect.ratio = 1)
lytic_wls_plot

ggsave(plot = lytic_wls_plot, filename = "./figures/lytic_WLS_bacteria_viruses.svg", dpi = 1000, width = 11, height = 7)

# WLS results ####
library(dplyr)
library(broom)

# Function to fit WLS for each predictor and extract model stats
wls_results <- long_lytic_df %>%
  group_by(Predictor) %>%
  group_modify(~ {
    # Fit Weighted Least Squares model
    wls_model <- lm(Corrected_VP_Lytic ~ Predictor_Value, data = .x, weights = .x$weights)
    
    # Extract model summary
    summary_wls <- summary(wls_model)
    
    # Store important results
    tibble(
      Intercept = coef(wls_model)[1],
      Estimate = coef(wls_model)[2],
      Std.Error = summary_wls$coefficients[2, 2],
      t_value = summary_wls$coefficients[2, 3],
      p_value = summary_wls$coefficients[2, 4],
      Adj_R2 = summary_wls$adj.r.squared,
      AIC = AIC(wls_model),
      n = nrow(.x)
    )
  })

# View results
print(wls_results)

write.csv(wls_results, "./results/stats/lytic_production_rate_WLS_bacteria_viruses.csv", row.names = F)


# 11. VBR Lytic production non-linear plot and results ####
# Compute weights for VBR using OLS residuals
vbr_lytic_df <- lytic_df %>%
  mutate(
    # Fit initial OLS model for VBR
    ols_model = list(lm(Corrected_VP_Lytic ~ VBR, data = .)),
    
    # Extract absolute residuals + small constant
    residuals_ols = abs(residuals(ols_model[[1]])) + 0.1, 
    
    # Compute weights as inverse of residuals
    weights = 1 / residuals_ols
  ) %>%
  dplyr::select(-ols_model, -residuals_ols)  # Remove intermediate columns

# Generate Weighted Regression Plot for VBR
vbr_wls_plot <- ggplot(vbr_lytic_df, aes(x = VBR, y = Corrected_VP_Lytic, shape = Season, fill = weights)) +
  # Weighted Robust Regression (Best Overall)
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = NA, se = T, 
              aes(weight = weights, group = 1), size = 1.2) +
  
  geom_point(size = 3,  stroke = 0.8, color = "black") +  # Circle points with black border
  scale_shape_manual(values = c(21, 24)) +  
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Weight") +  # Shared color scale
  
  # Weighted Regression (Best LM)
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#850000", se = F, 
              aes(weight = weights, group = 1), size = 1.2) +
  
  # Weighted Robust Regression (Best Overall)
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = "#0c1844", se = F, 
              aes(weight = weights, group = 1), size = 1.2) +
  
  # Legend for selected models
  annotate("text", x = 50, y = 9, 
           label = "Weighted", color = "#850000", hjust = 1, size = 5) +
  annotate("text", x = 50 , y = 8.5, 
           label = "Weighted Robust", color = "#0c1844", hjust = 1, size = 5) +
  
  # Labels and theme
  labs(x = "VBR", 
       y = expression("Lytic viral production rate " ~ (10^{5} ~ viruses ~ mL^{-1} ~ h^{-1})), 
       fill = "Weight") +
  
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        panel.border = element_rect(linewidth = 1.5, fill = NA),
        axis.line = element_line(size = 0), 
        axis.text = element_text(color = "black"),
        aspect.ratio = 1)

# Display the plot
vbr_wls_plot
# Save the plot
ggsave(plot = vbr_wls_plot, filename = "./figures/lytic_WLS_vbr.svg", dpi = 1000, width = 8, height = 7)





final_poly_table <- data.frame(
  Model = c("Unweighted", "Weighted", "Robust", "Weighted Robust"),
  Quadratic_Term = c(-4.250, -3.552, -3.349, -3.129),
  Std_Error = c(1.902, 1.190, 1.360, 1.179),
  CI_Lower = c(-8.552, -6.244, -6.016, -5.441),
  CI_Upper = c(0.052, -0.860, -0.683, -0.816),
  Adj_R2 = c(0.238, 0.474, NA, NA),
  AIC = c(54.03, 43.23, 54.46, 41.67),
  p_value = c(0.0523, 0.0153, NA, NA)  # NA for RLM models
)

# View table
final_poly_table
write.csv(final_poly_table, "./results/stats/lytic_production_rate_polynomial_VBR.csv", row.names = F)

