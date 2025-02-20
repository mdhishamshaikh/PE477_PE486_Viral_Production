# To explain lysogenic viral production rates #

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

lysogenic_df <- pe_df_7m %>%
  dplyr::filter(Corrected_VP_Lysogenic > 0)

# 1.0 Lysogenic viral production rate vs Total bacteria ####

bacteria_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ Total_Bacteria, data = lysogenic_df)

par(mfrow=c(2,2))
plot(bacteria_lysogenic_lm)


summary(bacteria_lysogenic_lm)

weights <- 1 / (abs(bacteria_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

bacteria_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ Total_Bacteria, 
                                 data = lysogenic_df, 
                                 weights = weights)

par(mfrow=c(2,2))
plot(bacteria_lysogenic_weighted_lm)
summary(bacteria_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(bacteria_lysogenic_weighted_lm)

AIC(bacteria_lysogenic_lm, bacteria_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(bacteria_lysogenic_lm)
summary_weighted <- summary(bacteria_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(bacteria_lysogenic_lm)[1], 2), 
                        " + ", round(coef(bacteria_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(bacteria_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(bacteria_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = Total_Bacteria, y = Corrected_VP_Lysogenic)) +
  geom_point(size = 4) +
  scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$Total_Bacteria), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$Total_Bacteria), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Total Bacteria", y = "Corrected VP Lysogenic") +
  theme_test()


# 2.0 Lysogenic viral production rate vs HNA ####

hna_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ HNA, data = lysogenic_df)

par(mfrow=c(2,2))
plot(hna_lysogenic_lm)

summary(hna_lysogenic_lm)

weights <- 1 / (abs(hna_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

hna_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ HNA, 
                            data = lysogenic_df, 
                            weights = weights)

par(mfrow=c(2,2))
plot(hna_lysogenic_weighted_lm)
summary(hna_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(hna_lysogenic_weighted_lm)

AIC(hna_lysogenic_lm, hna_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(hna_lysogenic_lm)
summary_weighted <- summary(hna_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(hna_lysogenic_lm)[1], 2), 
                        " + ", round(coef(hna_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(hna_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(hna_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = HNA, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$HNA), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$HNA), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "HNA bacteria", y = "Corrected VP Lysogenic") +
  theme_test()

# 3.0 Lysogenic viral production rate vs LNA ####

lna_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ LNA, data = lysogenic_df)

par(mfrow=c(2,2))
plot(lna_lysogenic_lm)

summary(lna_lysogenic_lm)

weights <- 1 / (abs(lna_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

lna_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ LNA, 
                            data = lysogenic_df, 
                            weights = weights)

par(mfrow=c(2,2))
plot(lna_lysogenic_weighted_lm)
summary(lna_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(lna_lysogenic_weighted_lm)

AIC(lna_lysogenic_lm, lna_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(lna_lysogenic_lm)
summary_weighted <- summary(lna_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(lna_lysogenic_lm)[1], 2), 
                        " + ", round(coef(lna_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(lna_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(lna_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = LNA, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$LNA), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$LNA), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "LNA bacteria", y = "Corrected VP Lysogenic") +
  theme_test()

# 4.0 Lysogenic viral production rate vs Total viruses ####

viruses_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ Total_Viruses, data = lysogenic_df )

par(mfrow=c(2,2))
plot(viruses_lysogenic_lm)

summary(viruses_lysogenic_lm)

weights <- 1 / (abs(viruses_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

viruses_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ Total_Viruses, 
                                data = lysogenic_df , 
                                weights = weights)

par(mfrow=c(2,2))
plot(viruses_lysogenic_weighted_lm)
summary(viruses_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(viruses_lysogenic_weighted_lm)

AIC(viruses_lysogenic_lm, viruses_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(viruses_lysogenic_lm)
summary_weighted <- summary(viruses_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(viruses_lysogenic_lm)[1], 2), 
                        " + ", round(coef(viruses_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(viruses_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(viruses_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = Total_Viruses, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$Total_Viruses), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$Total_Viruses), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Total Viruses", y = "Corrected VP Lysogenic") +
  theme_test()

# 5.0 Lysogenic viral production rate vs V1 ####

v1_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ V1, data = lysogenic_df)

par(mfrow=c(2,2))
plot(v1_lysogenic_lm)

summary(v1_lysogenic_lm)

weights <- 1 / (abs(v1_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

v1_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ V1, 
                           data = lysogenic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(v1_lysogenic_weighted_lm)
summary(v1_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v1_lysogenic_weighted_lm)

AIC(v1_lysogenic_lm, v1_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(v1_lysogenic_lm)
summary_weighted <- summary(v1_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v1_lysogenic_lm)[1], 2), 
                        " + ", round(coef(v1_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v1_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(v1_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = V1, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$V1), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$V1), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V1", y = "Corrected VP Lysogenic") +
  theme_test()

# 6.0 Lysogenic viral production rate vs V2 ####

v2_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ V2, data = lysogenic_df)

par(mfrow=c(2,2))
plot(v2_lysogenic_lm)

summary(v2_lysogenic_lm)

weights <- 1 / (abs(v2_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

v2_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ V2, 
                           data = lysogenic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(v2_lysogenic_weighted_lm)
summary(v2_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v2_lysogenic_weighted_lm)

AIC(v2_lysogenic_lm, v2_lysogenic_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(v2_lysogenic_lm)
summary_weighted <- summary(v2_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v2_lysogenic_lm)[1], 2), 
                        " + ", round(coef(v2_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v2_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(v2_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = V2, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$V2), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$V2), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V2", y = "Corrected VP Lysogenic") +
  theme_test()

# 7.0 Lysogenic viral production rate vs V3 ####

v3_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ V3, data = lysogenic_df)

par(mfrow=c(2,2))
plot(v3_lysogenic_lm)

summary(v3_lysogenic_lm)

weights <- 1 / (abs(v3_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

v3_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ V3, 
                           data = lysogenic_df, 
                           weights = weights)

par(mfrow=c(2,2))
plot(v3_lysogenic_weighted_lm)
summary(v3_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(v3_lysogenic_weighted_lm)

AIC(v3_lysogenic_lm, v3_lysogenic_weighted_lm)



# Get model summaries
summary_unweighted <- summary(v3_lysogenic_lm)
summary_weighted <- summary(v3_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(v3_lysogenic_lm)[1], 2), 
                        " + ", round(coef(v3_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(v3_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(v3_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = V3, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$V3), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$V3), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "V3", y = "Corrected VP Lysogenic") +
  theme_test()

# 8.0 Lysogenic viral production rate vs Cyanobacteria ####
cyano_df_filtered <- lysogenic_df %>% dplyr::filter(!Station_Number %in% c(5, 6))

cyano_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ Cyanobacteria, data = cyano_df_filtered)

par(mfrow=c(2,2))
plot(cyano_lysogenic_lm)

summary(cyano_lysogenic_lm)

weights <- 1 / (abs(cyano_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

cyano_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ Cyanobacteria, 
                              data = cyano_df_filtered, 
                              weights = weights)

par(mfrow=c(2,2))
plot(cyano_lysogenic_weighted_lm)
summary(cyano_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(cyano_lysogenic_weighted_lm)

AIC(cyano_lysogenic_lm, cyano_lysogenic_weighted_lm)


# Get model summaries
summary_unweighted <- summary(cyano_lysogenic_lm)
summary_weighted <- summary(cyano_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(cyano_lysogenic_lm)[1], 2), 
                        " + ", round(coef(cyano_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(cyano_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(cyano_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(cyano_df_filtered, aes(x = Cyanobacteria, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(cyano_df_filtered$Cyanobacteria), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(cyano_df_filtered$Cyanobacteria), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Cyanobacteria", y = "Corrected VP Lysogenic") +
  theme_test()

# 9.0 Lysogenic viral production with VBR ####



vbr_lysogenic_lm <- lm(Corrected_VP_Lysogenic ~ VBR, data = lysogenic_df)

par(mfrow=c(2,2))
plot(vbr_lysogenic_lm)

summary(vbr_lysogenic_lm)

weights <- 1 / (abs(vbr_lysogenic_lm$residuals) + 0.1)  # Avoid division by zero

vbr_lysogenic_weighted_lm <- lm(Corrected_VP_Lysogenic ~ VBR, 
                            data = lysogenic_df, 
                            weights = weights)

par(mfrow=c(2,2))
plot(vbr_lysogenic_weighted_lm)
summary(vbr_lysogenic_weighted_lm)
par(mfrow=c(1,1))
influencePlot(vbr_lysogenic_weighted_lm)

AIC(vbr_lysogenic_lm, vbr_lysogenic_weighted_lm)
# weighted regression is better.  

# Get model summaries
summary_unweighted <- summary(vbr_lysogenic_lm)
summary_weighted <- summary(vbr_lysogenic_weighted_lm)

# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(vbr_lysogenic_lm)[1], 2), 
                        " + ", round(coef(vbr_lysogenic_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(vbr_lysogenic_weighted_lm)[1], 2), 
                      " + ", round(coef(vbr_lysogenic_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(lysogenic_df, aes(x = VBR, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$VBR), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$VBR), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "VBR", y = "Corrected VP Lysogenic") +
  theme_test()


#POLYNOMIAL#
# Load required libraries
library(MASS)  # For robust regression
library(ggplot2)
library(car)   # For influence plots

# Fit an unweighted second-degree polynomial regression (LM)
vbr_lysogenic_poly <- lm(Corrected_VP_Lysogenic ~ poly(VBR, 2), data = lysogenic_df)

# Compute weights: Inverse of absolute residuals + small constant to avoid division by zero
weights <- 1 / (abs(vbr_lysogenic_poly$residuals) + 0.1)

# Fit weighted polynomial regression (LM)
vbr_lysogenic_weighted_poly <- lm(Corrected_VP_Lysogenic ~ poly(VBR, 2), 
                              data = lysogenic_df, 
                              weights = weights)

# Fit robust polynomial regression (RLM)
vbr_lysogenic_robust_poly <- rlm(Corrected_VP_Lysogenic ~ poly(VBR, 2), data = lysogenic_df)

# Compute new weights for robust regression
weights_robust <- 1 / (abs(resid(vbr_lysogenic_robust_poly)) + 0.1)

# Fit weighted robust polynomial regression (RLM with weights)
vbr_lysogenic_weighted_robust_poly <- rlm(Corrected_VP_Lysogenic ~ poly(VBR, 2), 
                                      data = lysogenic_df, 
                                      weights = weights_robust)

# Summaries
summary_unweighted <- summary(vbr_lysogenic_poly)
summary_weighted <- summary(vbr_lysogenic_weighted_poly)
summary_robust <- summary(vbr_lysogenic_robust_poly)
summary_weighted_robust <- summary(vbr_lysogenic_weighted_robust_poly)

# Compare AIC values
AIC(vbr_lysogenic_poly, vbr_lysogenic_weighted_poly, vbr_lysogenic_robust_poly, vbr_lysogenic_weighted_robust_poly)

# Diagnostic Plots
par(mfrow=c(2,2))
plot(vbr_lysogenic_poly)
plot(vbr_lysogenic_weighted_poly)
plot(vbr_lysogenic_robust_poly)
plot(vbr_lysogenic_weighted_robust_poly)

# Influence plots
par(mfrow=c(1,1))
influencePlot(vbr_lysogenic_weighted_poly)
influencePlot(vbr_lysogenic_weighted_robust_poly)

# Extract model equations and statistics
eq_unweighted <- paste0("y = ", round(coef(vbr_lysogenic_poly)[1], 2), 
                        " + ", round(coef(vbr_lysogenic_poly)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(vbr_lysogenic_weighted_poly)[1], 2), 
                      " + ", round(coef(vbr_lysogenic_weighted_poly)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

eq_robust <- paste0("y = ", round(coef(vbr_lysogenic_robust_poly)[1], 2), 
                    " + ", round(coef(vbr_lysogenic_robust_poly)[2], 2), "*x\n",
                    "Adj R² = NA (RLM)", 
                    ", p = NA (RLM)")

eq_weighted_robust <- paste0("y = ", round(coef(vbr_lysogenic_weighted_robust_poly)[1], 2), 
                             " + ", round(coef(vbr_lysogenic_weighted_robust_poly)[2], 2), "*x\n",
                             "Adj R² = NA (RLM)", 
                             ", p = NA (RLM)")

# Create scatterplots with polynomial regression lines
ggplot(lysogenic_df, aes(x = VBR, y = Corrected_VP_Lysogenic)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#0c1844", se = FALSE) +
  annotate("text", x = min(lysogenic_df$VBR), y = 12, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#850000", se = FALSE, 
              aes(weight = weights)) +
  annotate("text", x = min(lysogenic_df$VBR), y = 10, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = "#008500", se = FALSE) +
  annotate("text", x =35, y = 12, 
           label = eq_robust, hjust = 0, size = 5, color = "#008500") +
  geom_smooth(method = "rlm", formula = y ~ poly(x, 2), color = "#856300", se = FALSE, 
              aes(weight = weights_robust)) +
  annotate("text", x = 35, y = 10, 
           label = eq_weighted_robust, hjust = 0, size = 5, color = "#856300") +
  labs(#title = "Polynomial Regression: Unweighted vs. Weighted vs. Robust",
    x = "VBR", y = "Corrected VP Lysogenic") +
  theme_test()


confint(vbr_lysogenic_poly)
confint(vbr_lysogenic_weighted_poly)
confint(vbr_lysogenic_robust_poly)  
confint(vbr_lysogenic_weighted_robust_poly)  


