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
                decay_rate_linear  = decay_rate_linear /1e+5,
                Corrected_VP_Lytic  = Corrected_VP_Lytic/1e+5,
                Corrected_VP_Lysogenic = Corrected_VP_Lysogenic/1e+5
  )

pe_df_7m <- pe_df %>%
  dplyr::filter(Depth == 7)
decay_df <- pe_df_7m %>%
  dplyr::filter(decay_rate_linear > 0)

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



# 3.0 Viral decay rate vs Max_Depth####

max_depth_decay_lm <- lm(decay_rate_linear ~ Max_Depth, data = decay_df)

par(mfrow=c(2,2))
plot(max_depth_decay_lm)
par(mfrow=c(1,1))
influencePlot(max_depth_decay_lm)
summary(max_depth_decay_lm)

weights <- 1 / (abs(max_depth_decay_lm$residuals) + 0.1)  # Avoid division by zero

max_depth_decay_weighted_lm <- lm(decay_rate_linear ~ Max_Depth, 
                                          data = decay_df, 
                                          weights = weights)

par(mfrow=c(2,2))
plot(max_depth_decay_weighted_lm)
summary(max_depth_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(max_depth_decay_weighted_lm)

AIC(max_depth_decay_lm, max_depth_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(max_depth_decay_lm)


summary_weighted <- summary(max_depth_decay_weighted_lm)

# Neither assuptions were met.
# Trying polynomial second degree





# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(max_depth_decay_lm)[1], 2), 
                              " + ", round(coef(max_depth_decay_lm)[2], 2), "*x\n",
                              "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                              ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(max_depth_decay_weighted_lm)[1], 2), 
                            " + ", round(coef(max_depth_decay_weighted_lm)[2], 2), "*x\n",
                            "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                            ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(decay_df, aes(x = Max_Depth, y = decay_rate_linear, shape = Season)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = F,
              aes(group = 1)) +
  annotate("text", x = min(decay_df$Max_Depth), y = 240, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights, group = 1)) +
  annotate("text", x = min(decay_df$Max_Depth), y = 200, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Max_Depth", y = "Viral decay rate") +
  theme_test()





# 4.0 Viral decay rate vs Total_Viruses####

viruses_decay_lm <- lm(decay_rate_linear ~ Total_Viruses, data = decay_df)

par(mfrow=c(2,2))
plot(viruses_decay_lm)
par(mfrow=c(1,1))
influencePlot(viruses_decay_lm)
summary(viruses_decay_lm)

weights <- 1 / (abs(viruses_decay_lm$residuals) + 0.1)  # Avoid division by zero

viruses_decay_weighted_lm <- lm(decay_rate_linear ~ Total_Viruses, 
                                  data = decay_df, 
                                  weights = weights)

par(mfrow=c(2,2))
plot(viruses_decay_weighted_lm)
summary(viruses_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(viruses_decay_weighted_lm)

AIC(viruses_decay_lm, viruses_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(viruses_decay_lm)


summary_weighted <- summary(viruses_decay_weighted_lm)

# Neither assuptions were met.
# Trying polynomial second degree





# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(viruses_decay_lm)[1], 2), 
                        " + ", round(coef(viruses_decay_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(viruses_decay_weighted_lm)[1], 2), 
                      " + ", round(coef(viruses_decay_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(decay_df, aes(x = Total_Viruses, y = decay_rate_linear, shape = Season)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = F,
              aes(group = 1)) +
  annotate("text", x = min(decay_df$Total_Viruses), y = 280, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights, group = 1)) +
  annotate("text", x = min(decay_df$Total_Viruses), y = 240, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Total_Viruses", y = "Viral decay rate") +
  theme_test()


# 4.0 Viral percent decay vs Total_Viruses####

viruses_percent_decay_lm <- lm(percent_decay_day_linear ~ Total_Viruses, data = decay_df)

par(mfrow=c(2,2))
plot(viruses_percent_decay_lm)
par(mfrow=c(1,1))
influencePlot(viruses_percent_decay_lm)
summary(viruses_percent_decay_lm)

weights <- 1 / (abs(viruses_percent_decay_lm$residuals) + 0.1)  # Avoid division by zero

viruses_percent_decay_weighted_lm <- lm(percent_decay_day_linear ~ Total_Viruses, 
                                data = decay_df, 
                                weights = weights)

par(mfrow=c(2,2))
plot(viruses_percent_decay_weighted_lm)
summary(viruses_percent_decay_weighted_lm)
par(mfrow=c(1,1))
influencePlot(viruses_percent_decay_weighted_lm)

AIC(viruses_percent_decay_lm, viruses_percent_decay_weighted_lm)
# weighted regression is better.  


# Get model summaries
summary_unweighted <- summary(viruses_percent_decay_lm)


summary_weighted <- summary(viruses_percent_decay_weighted_lm)

# Neither assuptions were met.
# Trying polynomial second degree





# Extract model statistics
eq_unweighted <- paste0("y = ", round(coef(viruses_percent_decay_lm)[1], 2), 
                        " + ", round(coef(viruses_percent_decay_lm)[2], 2), "*x\n",
                        "Adj R² = ", round(summary_unweighted$adj.r.squared, 3), 
                        ", p = ", round(summary_unweighted$coefficients[2,4], 3))

eq_weighted <- paste0("y = ", round(coef(viruses_percent_decay_weighted_lm)[1], 2), 
                      " + ", round(coef(viruses_percent_decay_weighted_lm)[2], 2), "*x\n",
                      "Adj R² = ", round(summary_weighted$adj.r.squared, 3), 
                      ", p = ", round(summary_weighted$coefficients[2,4], 3))

# Create scatterplots
ggplot(decay_df, aes(x = Total_Viruses, y = percent_decay_day_linear, shape = Season)) +
  geom_point(size = 4) +
  #scale_shape_manual(values = custom_shape_palette_cruise) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = F,
              aes(group = 1)) +
  annotate("text", x = min(decay_df$Total_Viruses), y = 50, 
           label = eq_unweighted, hjust = 0, size = 5, color = "#0c1844") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = FALSE, 
              aes(weight = weights, group = 1)) +
  annotate("text", x = min(decay_df$Total_Viruses), y = 40, 
           label = eq_weighted, hjust = 0, size = 5, color = "#850000") +
  labs(#title = "Unweighted Regression", 
    x = "Total_Viruses", y = "Viral decay percent") +
  theme_test()














# 4.0 Plots for Viral decay, turbiidty and temperature #####

# Here I will make a combined (both seasons) and per season weighted linear regression

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Create Subsets
decay_complete <- decay_df %>% mutate(Dataset = "Complete")
decay_autumn <- decay_df %>% filter(Season == "Autumn (Sept 2020)") %>% mutate(Dataset = "Autumn")
decay_spring <- decay_df %>% filter(Season == "Spring (Apr 2021)") %>% mutate(Dataset = "Spring")

# Step 2: Function to Compute Weights for Each Predictor in a Given Dataset
compute_weights <- function(df, predictor) {
  df %>%
    mutate(
      ols_model = list(lm(decay_rate_linear ~ .data[[predictor]], data = cur_data())),
      residuals_ols = abs(residuals(ols_model[[1]])) + 0.1,  # Avoid division by zero
      weights = 1 / residuals_ols
    ) %>%
    dplyr::select(-ols_model, -residuals_ols)
}

# Apply function for both predictors in each dataset
complete_temp <- compute_weights(decay_complete, "Temperature")
complete_turb <- compute_weights(decay_complete, "Turbidity")

autumn_temp <- compute_weights(decay_autumn, "Temperature")
autumn_turb <- compute_weights(decay_autumn, "Turbidity")

spring_temp <- compute_weights(decay_spring, "Temperature")
spring_turb <- compute_weights(decay_spring, "Turbidity")

# Step 3: Combine Data into a Unified DataFrame
combine_long_df <- function(df, predictor) {
  df %>%
    pivot_longer(cols = c("Temperature", "Turbidity"),
                 names_to = "Predictor",
                 values_to = "Predictor_Value") %>%
    filter(Predictor == predictor) %>%
    mutate(Predictor = factor(Predictor, levels = c("Temperature", "Turbidity")),
           Dataset = factor(Dataset, levels = c("Complete", "Autumn", "Spring")))
}

long_df <- bind_rows(
  combine_long_df(complete_temp, "Temperature"),
  combine_long_df(complete_turb, "Turbidity"),
  combine_long_df(autumn_temp, "Temperature"),
  combine_long_df(autumn_turb, "Turbidity"),
  combine_long_df(spring_temp, "Temperature"),
  combine_long_df(spring_turb, "Turbidity")
)

# Step 4: Faceted WLS Plot
decay_wls_turbidity_plot <- ggplot(long_df %>% dplyr::filter(Predictor == "Turbidity"), aes(x = Predictor_Value, y = decay_rate_linear, shape = Season, fill = weights)) +
  geom_smooth(method = "lm", formula = y ~ x, color = NA, se = T, size = 1, 
              aes(weight = weights, group = 1)) +
  geom_point(size = 3, stroke = 0.8, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Weight") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = F, size = 1, 
              aes(weight = weights, group = 1)) +
  facet_wrap(Predictor ~ Dataset, scales = "free_x", strip.position = "bottom",
             labeller = as_labeller(variable_labels, label_parsed)) +
  labs(x = NULL, 
       y = expression(atop("Viral decay rate", (10^{5} ~ viruses ~ mL^{-1} ~ h^{-1}))),
       #title = "Weighted Least Squares Regression (WLS) for Decay Rate",
       fill = "Weight") +
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
decay_wls_turbidity_plot
# Step 5: Save the Plot
ggsave(plot = decay_wls_turbidity_plot, filename = "./figures/decay_WLS_turbidity.svg", dpi = 1000, width = 11, height = 4)

# Step 4: Faceted WLS Plot
decay_wls_temperature_plot <- ggplot(long_df %>% dplyr::filter(Predictor == "Temperature"), aes(x = Predictor_Value, y = decay_rate_linear, shape = Season, fill = weights)) +
  geom_smooth(method = "lm", formula = y ~ x, color = NA, se = T, size = 1, 
              aes(weight = weights, group = 1)) +
  geom_point(size = 3, stroke = 0.8, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Weight") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = F, size = 1, 
              aes(weight = weights, group = 1)) +
  facet_wrap(Predictor ~ Dataset, scales = "free_x", strip.position = "bottom",
             labeller = as_labeller(variable_labels, label_parsed)) +
  labs(x = NULL, 
       y = expression(atop("Viral decay rate", (10^{5} ~ viruses ~ mL^{-1} ~ h^{-1}))),
       #title = "Weighted Least Squares Regression (WLS) for Decay Rate",
       fill = "Weight") +
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
decay_wls_temperature_plot
# Step 5: Save the Plot
ggsave(plot = decay_wls_temperature_plot, filename = "./figures/decay_WLS_temperature.svg", dpi = 1000, width = 11, height = 4)

# Function to fit WLS for each Predictor and Dataset, and extract model stats
decay_temp_turb_wls_results <- long_df %>%
  group_by(Predictor, Dataset) %>%
  group_modify(~ {
    # Fit Weighted Least Squares model
    wls_model <- lm(decay_rate_linear ~ Predictor_Value, data = .x, weights = .x$weights)
    
    # Extract model summary
    summary_wls <- summary(wls_model)
    
    # Store important results
    tibble(
      Intercept = coef(wls_model)[1],
      Estimate = coef(wls_model)[2],
      Std_Error = summary_wls$coefficients[2, 2],
      t_value = summary_wls$coefficients[2, 3],
      p_value = summary_wls$coefficients[2, 4],
      Adj_R2 = summary_wls$adj.r.squared,
      AIC = AIC(wls_model),
      n = nrow(.x)
    )
  })

# Print WLS results
print(decay_temp_turb_wls_results)
write.csv(decay_temp_turb_wls_results, "./results/stats/decay_rate_WLS_temperature_turbidity.csv", row.names = F)


par(mfrow=c(2,2))
plot(lm(decay_rate_linear ~ Temperature, data = complete_temp, weights = complete_temp$weights))
plot(lm(decay_rate_linear ~ Temperature, data = autumn_temp, weights = autumn_temp$weights))
plot(lm(decay_rate_linear ~ Temperature, data = spring_temp, weights = spring_temp$weights))
plot(lm(decay_rate_linear ~ Turbidity, data = complete_turb, weights = complete_turb$weights))
plot(lm(decay_rate_linear ~ Turbidity, data = autumn_turb, weights = autumn_turb$weights))
plot(lm(decay_rate_linear ~ Turbidity, data = spring_turb, weights = spring_turb$weights))


# 5.0 pLOTS FOR Viral decay and max depth ####
# Linear
# Weighted linear
# Log linear
# Plynomial - second degree

# Load necessary libraries
library(ggplot2)
library(car)  # For influencePlot

# Step 1: Fit Models
## 1. Linear Model
max_depth_decay_lm <- lm(decay_rate_linear ~ Max_Depth, data = decay_df)

## 2. Weighted Linear Model
weights <- 1 / (abs(residuals(max_depth_decay_lm)) + 0.1)  # Avoid division by zero
max_depth_decay_weighted_lm <- lm(decay_rate_linear ~ Max_Depth, data = decay_df, weights = weights)

## 3. Log-Linear Model
max_depth_decay_log_lm <- lm(log(decay_rate_linear) ~ Max_Depth, data = decay_df)

## 4. Polynomial Model (Quadratic)
max_depth_decay_poly_lm <- lm(decay_rate_linear ~ poly(Max_Depth, 2), data = decay_df)

# Step 2: Model Diagnostics (Built-in `plot()` function)
par(mfrow=c(2,2))
plot(max_depth_decay_lm, main="Linear Model")

par(mfrow=c(2,2))
plot(max_depth_decay_weighted_lm, main="Weighted Linear Model")

par(mfrow=c(2,2))
plot(max_depth_decay_log_lm, main="Log-Linear Model")

par(mfrow=c(2,2))
plot(max_depth_decay_poly_lm, main="Polynomial Model")

# Reset plotting layout
par(mfrow=c(1,1))

# Step 3: Compare Models using AIC
aic_results <- AIC(
  max_depth_decay_lm,
  max_depth_decay_weighted_lm,
  max_depth_decay_log_lm,
  max_depth_decay_poly_lm
)
print(aic_results)

# Step 4: Create Comparison Plot
ggplot(decay_df, aes(x = Max_Depth, y = decay_rate_linear, shape = Season)) +
  geom_point(size = 4) +
  
  # Linear Model (Unweighted)
  geom_smooth(method = "lm", formula = y ~ x, color = "#0c1844", se = F, aes(group = 1)) +
  
  # Weighted Linear Model
  geom_smooth(method = "lm", formula = y ~ x, color = "#850000", se = F, aes(weight = weights, group = 1)) +
  
  # Log-Linear Model
  geom_smooth(method = "lm", formula = y ~ log(x), color = "#006400", se = F, aes(group = 1)) +
  
  # Polynomial Model
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#FF8C00", se = F, aes(group = 1)) +
  
  labs(x = "Max Depth", y = "Viral Decay Rate") +
  theme_test()


# Step 2: Function to Extract Model Statistics
# Load necessary libraries
library(dplyr)
library(tibble)

# Load necessary libraries
library(dplyr)
library(tibble)

# Function to Extract Model Statistics with Confidence Intervals
extract_model_stats <- function(model, model_name, poly = FALSE) {
  summary_model <- summary(model)
  model_coefs <- coef(summary_model)
  conf_int <- confint(model, level = 0.95)  # Compute 95% confidence intervals
  
  if (poly) {
    # Extract both linear and quadratic terms separately
    tibble(
      Model = model_name,
      Intercept_Estimate = round(model_coefs[1, 1], 2),
      Intercept_SE = round(model_coefs[1, 2], 2),
      Intercept_Lower_CI = round(conf_int[1, 1], 2),
      Intercept_Upper_CI = round(conf_int[1, 2], 2),
      
      Linear_Term_Estimate = round(model_coefs[2, 1], 2),
      Linear_Term_SE = round(model_coefs[2, 2], 2),
      Linear_Term_Lower_CI = round(conf_int[2, 1], 2),
      Linear_Term_Upper_CI = round(conf_int[2, 2], 2),
      
      Quadratic_Term_Estimate = round(model_coefs[3, 1], 2),
      Quadratic_Term_SE = round(model_coefs[3, 2], 2),
      Quadratic_Term_Lower_CI = round(conf_int[3, 1], 2),
      Quadratic_Term_Upper_CI = round(conf_int[3, 2], 2),
      
      t_Linear = round(model_coefs[2, 3], 2),
      p_Linear = round(model_coefs[2, 4], 5),
      t_Quadratic = round(model_coefs[3, 3], 2),
      p_Quadratic = round(model_coefs[3, 4], 5),
      
      Adj_R2 = round(summary_model$adj.r.squared, 3),
      AIC = round(AIC(model), 2),
      Residual_SE = round(summary_model$sigma, 3)
    )
  } else {
    # Extract only the linear term separately
    tibble(
      Model = model_name,
      Intercept_Estimate = round(model_coefs[1, 1], 2),
      Intercept_SE = round(model_coefs[1, 2], 2),
      Intercept_Lower_CI = round(conf_int[1, 1], 2),
      Intercept_Upper_CI = round(conf_int[1, 2], 2),
      
      Linear_Term_Estimate = round(model_coefs[2, 1], 2),
      Linear_Term_SE = round(model_coefs[2, 2], 2),
      Linear_Term_Lower_CI = round(conf_int[2, 1], 2),
      Linear_Term_Upper_CI = round(conf_int[2, 2], 2),
      
      Quadratic_Term_Estimate = NA,
      Quadratic_Term_SE = NA,
      Quadratic_Term_Lower_CI = NA,
      Quadratic_Term_Upper_CI = NA,
      
      t_Linear = round(model_coefs[2, 3], 2),
      p_Linear = round(model_coefs[2, 4], 5),
      t_Quadratic = NA,
      p_Quadratic = NA,
      
      Adj_R2 = round(summary_model$adj.r.squared, 3),
      AIC = round(AIC(model), 2),
      Residual_SE = round(summary_model$sigma, 3)
    )
  }
}

# Extract Model Summaries with Confidence Intervals
deacy_max_depth_model_results <- bind_rows(
  extract_model_stats(max_depth_decay_lm, "Unweighted Linear"),
  extract_model_stats(max_depth_decay_weighted_lm, "Weighted Linear"),
  extract_model_stats(max_depth_decay_log_lm, "Log-Linear"),
  extract_model_stats(max_depth_decay_poly_lm, "Polynomial (Quadratic)", poly = TRUE)
)

# Print Model Results
print(deacy_max_depth_model_results)

write.csv(deacy_max_depth_model_results, "./results/stats/decay_rate_poly_max_depth.csv", row.names = F)


# Compute 95% Confidence Intervals for Polynomial Model
confint(max_depth_decay_poly_lm, level = 0.95)

confint(max_depth_decay_lm, level = 0.95)
# Load necessary libraries
library(ggplot2)
library(viridis)

# Polynomial Regression Plot with Styled Theme
decay_poly_plot <- ggplot(decay_df, aes(x = Max_Depth, y = decay_rate_linear, shape = Season)) +
  
  # Confidence Interval (SE)
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = NA, se = TRUE, size = 1, aes(group = 1)) +
  
  # Data Points
  geom_point(size = 3, stroke = 0.8, color = "black") +
  
  
  # Polynomial Model Fit (Black Line)
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", se = FALSE, size = 1, aes(group = 1)) +
  
  # Labels & Theme
  labs(x = "Maximum depth (m)", 
       y = expression(atop("Viral decay rate", (10^{5} ~ viruses ~ mL^{-1} ~ h^{-1})))
       ) +
  
  # Custom Theme (Matches Your WLS Temperature Plot)
  theme_test(base_size = 15) +
  theme(legend.position = "right",
        panel.spacing.y = unit(0.5, "cm"),
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

# Print the plot
decay_poly_plot
ggsave(plot = decay_poly_plot, filename = "./figures/decay_poly_max_depth.svg", dpi = 1000, width = 11, height = 4)

