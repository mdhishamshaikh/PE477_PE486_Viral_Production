# Setting up 

source("./scripts/0_source.R")


# Loading CTD profiles
ctd_profiles <- read.csv("./results/ctd_profiles/ctd_profiles.csv")


# Calculating potential temperature using OCE R package
ctd_profiles <- ctd_profiles %>%
  mutate(Potential_tempertaure =oce::swTheta(salinity = Salinity, temperature = Temperature, pressure = Pressure))


plot( ctd_profiles$Temperature, ctd_profiles$Potential_tempertaure)
# Practically the same. I will use the temperature for now


# using Davidatlarge's ggTS function
# https://github.com/Davidatlarge/ggTS/blob/master/ggTS_DK.R


# Function to create a ggTS plot for a given subset ####
create_TS_plot <- function(subset_data, location_station) {
  ggTS(
    sal = subset_data$Salinity,
    pot.temp = subset_data$Temperature,
    reference.p = 0,
    col.par = subset_data$Depth,
    col.name = "depth [m]"
  ) +
    labs(title = paste(location_station)) +
    ylab("Temperature [°C]") +
    xlab("Salinity [psu]") +
    scale_color_viridis_c(option = "inferno", direction = -1, name = "Depth [m]", guide = guide_colorbar(reverse = TRUE)) +
    theme_bw() +
    theme(
      legend.position = "right"  # Ensure the legend is displayed on the right
    )
}

# Split data by Location_Station_Number
ts_plots <- ctd_profiles %>%
  split(.$Location_Station_Number) %>%
  imap(~ create_TS_plot(.x, .y)) 

# Combine all plots using cowplot
combined_plot <- plot_grid(plotlist = ts_plots, ncol = 4)

# Save combined plot
ggsave("./figures/combined_TS_plots.svg", plot = combined_plot, width = 14, height = 10, dpi = 800)
