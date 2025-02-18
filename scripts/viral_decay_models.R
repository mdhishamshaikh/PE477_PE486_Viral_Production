# To explain viral decay rates #

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

pe477_decay_df <- pe_df_7m %>%
  dplyr::filter(Location == "PE477")
pe486_decay_df <- pe_df_7m %>%
  dplyr::filter(Location == "PE486")


# 1.0 Viral decay rate vs Turbidity####
# 1.1 pe477
pe477_turbidity_decay_lm <- lm(decay_rate_linear ~ Turbidity, data = pe477_decay_df)

par(mfrow=c(2,2))
plot(pe477_turbidity_decay_lm)

summary(pe477_turbidity_decay_lm)

pe477_weights <- 1 / (abs(pe477_turbidity_decay_lm$residuals) + 0.1)  # Avoid division by zero

pe477_turbidity_decay_weighted_lm <- lm(decay_rate_linear ~ Turbidity, 
                                 data = pe477_decay_df, 
                                 weights = pe477_weights)

par(mfrow=c(2,2))
plot(pe477_turbidity_decay_weighted_lm)
summary(pe477_turbidity_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(pe477_turbidity_decay_weighted_lm)

AIC(pe477_turbidity_decay_lm, pe477_turbidity_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
pe477_summary_unweighted <- summary(pe477_turbidity_decay_lm)
pe477_summary_weighted <- summary(pe477_turbidity_decay_weighted_lm)

# Extract model statistics
pe477_eq_unweighted <- paste0("y = ", round(coef(pe477_turbidity_decay_lm)[1], 2), 
                        " + ", round(coef(pe477_turbidity_decay_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(pe477_summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(pe477_summary_unweighted$coefficients[2,4], 3))

pe477_eq_weighted <- paste0("y = ", round(coef(pe477_turbidity_decay_weighted_lm)[1], 2), 
                      " + ", round(coef(pe477_turbidity_decay_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(pe477_summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(pe477_summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(pe477_decay_df, aes(x = Turbidity, y = decay_rate_linear)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(pe477_decay_df$Turbidity), y = 240, 
           label = pe477_eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = pe477_weights)) +
  annotate("text", x = min(pe477_decay_df$Turbidity), y = 200, 
           label = pe477_eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Turbidity", y = "Viral decay rate") +
  theme_test()

# 1.2 pe486
pe486_turbidity_decay_lm <- lm(decay_rate_linear ~ Turbidity, data = pe486_decay_df)

par(mfrow=c(2,2))
plot(pe486_turbidity_decay_lm)

summary(pe486_turbidity_decay_lm)

pe486_weights <- 1 / (abs(pe486_turbidity_decay_lm$residuals) + 0.1)  # Avoid division by zero

pe486_turbidity_decay_weighted_lm <- lm(decay_rate_linear ~ Turbidity, 
                                        data = pe486_decay_df, 
                                        weights = pe486_weights)

par(mfrow=c(2,2))
plot(pe486_turbidity_decay_weighted_lm)
summary(pe486_turbidity_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(pe486_turbidity_decay_weighted_lm)

AIC(pe486_turbidity_decay_lm, pe486_turbidity_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
pe486_summary_unweighted <- summary(pe486_turbidity_decay_lm)
pe486_summary_weighted <- summary(pe486_turbidity_decay_weighted_lm)

# Extract model statistics
pe486_eq_unweighted <- paste0("y = ", round(coef(pe486_turbidity_decay_lm)[1], 2), 
                              " + ", round(coef(pe486_turbidity_decay_lm)[2], 2), "*x\n",
                              "Adj R² = ", round(pe486_summary_unweighted$adj.r.squared, 3), 
                              ", p = ", round(pe486_summary_unweighted$coefficients[2,4], 3))

pe486_eq_weighted <- paste0("y = ", round(coef(pe486_turbidity_decay_weighted_lm)[1], 2), 
                            " + ", round(coef(pe486_turbidity_decay_weighted_lm)[2], 2), "*x\n",
                            "Adj R² = ", round(pe486_summary_weighted$adj.r.squared, 3), 
                            ", p = ", round(pe486_summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(pe486_decay_df, aes(x = Turbidity, y = decay_rate_linear)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(pe486_decay_df$Turbidity), y = 280, 
           label = pe486_eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = pe486_weights)) +
  annotate("text", x = min(pe486_decay_df$Turbidity), y = 240, 
           label = pe486_eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Turbidity", y = "Viral decay rate") +
  theme_test()


# 2.0 Viral decay rate vs Temperature####
# 2.1 pe477
pe477_temperature_decay_lm <- lm(decay_rate_linear ~ Temperature, data = pe477_decay_df)

par(mfrow=c(2,2))
plot(pe477_temperature_decay_lm)

summary(pe477_temperature_decay_lm)

pe477_weights <- 1 / (abs(pe477_temperature_decay_lm$residuals) + 0.1)  # Avoid division by zero

pe477_temperature_decay_weighted_lm <- lm(decay_rate_linear ~ Temperature, 
                                        data = pe477_decay_df, 
                                        weights = pe477_weights)

par(mfrow=c(2,2))
plot(pe477_temperature_decay_weighted_lm)
summary(pe477_temperature_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(pe477_temperature_decay_weighted_lm)

AIC(pe477_temperature_decay_lm, pe477_temperature_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
pe477_summary_unweighted <- summary(pe477_temperature_decay_lm)
pe477_summary_weighted <- summary(pe477_temperature_decay_weighted_lm)

# Extract model statistics
pe477_eq_unweighted <- paste0("y = ", round(coef(pe477_temperature_decay_lm)[1], 2), 
                              " + ", round(coef(pe477_temperature_decay_lm)[2], 2), "*x\n",
                              "Adj R² = ", round(pe477_summary_unweighted$adj.r.squared, 3), 
                              ", p = ", round(pe477_summary_unweighted$coefficients[2,4], 3))

pe477_eq_weighted <- paste0("y = ", round(coef(pe477_temperature_decay_weighted_lm)[1], 2), 
                            " + ", round(coef(pe477_temperature_decay_weighted_lm)[2], 2), "*x\n",
                            "Adj R² = ", round(pe477_summary_weighted$adj.r.squared, 3), 
                            ", p = ", round(pe477_summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(pe477_decay_df, aes(x = Temperature, y = decay_rate_linear)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(pe477_decay_df$Temperature), y = 240, 
           label = pe477_eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = pe477_weights)) +
  annotate("text", x = min(pe477_decay_df$Temperature), y = 200, 
           label = pe477_eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Temperature", y = "Viral decay rate") +
  theme_test()

# 2.2 pe486
pe486_temperature_decay_lm <- lm(decay_rate_linear ~ Temperature, data = pe486_decay_df)

par(mfrow=c(2,2))
plot(pe486_temperature_decay_lm)

summary(pe486_temperature_decay_lm)

pe486_weights <- 1 / (abs(pe486_temperature_decay_lm$residuals) + 0.1)  # Avoid division by zero

pe486_temperature_decay_weighted_lm <- lm(decay_rate_linear ~ Temperature, 
                                        data = pe486_decay_df, 
                                        weights = pe486_weights)

par(mfrow=c(2,2))
plot(pe486_temperature_decay_weighted_lm)
summary(pe486_temperature_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(pe486_temperature_decay_weighted_lm)

AIC(pe486_temperature_decay_lm, pe486_temperature_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
pe486_summary_unweighted <- summary(pe486_temperature_decay_lm)
pe486_summary_weighted <- summary(pe486_temperature_decay_weighted_lm)

# Extract model statistics
pe486_eq_unweighted <- paste0("y = ", round(coef(pe486_temperature_decay_lm)[1], 2), 
                              " + ", round(coef(pe486_temperature_decay_lm)[2], 2), "*x\n",
                              "Adj R² = ", round(pe486_summary_unweighted$adj.r.squared, 3), 
                              ", p = ", round(pe486_summary_unweighted$coefficients[2,4], 3))

pe486_eq_weighted <- paste0("y = ", round(coef(pe486_temperature_decay_weighted_lm)[1], 2), 
                            " + ", round(coef(pe486_temperature_decay_weighted_lm)[2], 2), "*x\n",
                            "Adj R² = ", round(pe486_summary_weighted$adj.r.squared, 3), 
                            ", p = ", round(pe486_summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(pe486_decay_df, aes(x = Temperature, y = decay_rate_linear)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(pe486_decay_df$Temperature), y = 280, 
           label = pe486_eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = pe486_weights)) +
  annotate("text", x = min(pe486_decay_df$Temperature), y = 240, 
           label = pe486_eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Temperature", y = "Viral decay rate") +
  theme_test()



# 3.0 Viral decay rate vs Temperature####
# 2.1 pe477
pe477_temperature_decay_lm <- lm(decay_rate_linear ~ Temperature, data = pe477_decay_df)

par(mfrow=c(2,2))
plot(pe477_temperature_decay_lm)

summary(pe477_temperature_decay_lm)

pe477_weights <- 1 / (abs(pe477_temperature_decay_lm$residuals) + 0.1)  # Avoid division by zero

pe477_temperature_decay_weighted_lm <- lm(decay_rate_linear ~ Temperature, 
                                          data = pe477_decay_df, 
                                          weights = pe477_weights)

par(mfrow=c(2,2))
plot(pe477_temperature_decay_weighted_lm)
summary(pe477_temperature_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(pe477_temperature_decay_weighted_lm)

AIC(pe477_temperature_decay_lm, pe477_temperature_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
pe477_summary_unweighted <- summary(pe477_temperature_decay_lm)
pe477_summary_weighted <- summary(pe477_temperature_decay_weighted_lm)

# Extract model statistics
pe477_eq_unweighted <- paste0("y = ", round(coef(pe477_temperature_decay_lm)[1], 2), 
                              " + ", round(coef(pe477_temperature_decay_lm)[2], 2), "*x\n",
                              "Adj R² = ", round(pe477_summary_unweighted$adj.r.squared, 3), 
                              ", p = ", round(pe477_summary_unweighted$coefficients[2,4], 3))

pe477_eq_weighted <- paste0("y = ", round(coef(pe477_temperature_decay_weighted_lm)[1], 2), 
                            " + ", round(coef(pe477_temperature_decay_weighted_lm)[2], 2), "*x\n",
                            "Adj R² = ", round(pe477_summary_weighted$adj.r.squared, 3), 
                            ", p = ", round(pe477_summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(pe477_decay_df, aes(x = Temperature, y = decay_rate_linear)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = FALSE) +
  annotate("text", x = min(pe477_decay_df$Temperature), y = 240, 
           label = pe477_eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = pe477_weights)) +
  annotate("text", x = min(pe477_decay_df$Temperature), y = 200, 
           label = pe477_eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Temperature", y = "Viral decay rate") +
  theme_test()


