library(tidyverse)
library(leaflet)
library(ggsci)
library(ggmap)
library(sf)
library(cowplot)
library(patchwork)

pe_df<- read.csv("./PE_Cruises/results/PE_Cruises_viral_production/vp_abundance_nutrients_ts.csv")
pe_df <- pe_df %>%
  mutate(Sample_Index = factor(paste0(Location, "-", Station_Number)))

#Let's define some categories.
nutrients<- c("TON", "Nitrite", #'Nitrate', #remove nitrate as it is claculated using TON-Nitrite
              "Phosphate", "Silicate")
ts<- c("Temperature",
       "Salinity")
viruses<- c("Total_Viruses",
            "V1", "V2", "V3")
bacteria<- c("Total_Bacteria", "HNA", "LNA")
common<- c("Location", "Station_Number", 
           "Depth", "Sample_Index")


#Nutrients ####

nuts_df <- pe_df %>%
  select(all_of(c(common, nutrients))) %>%
  pivot_longer(cols = nutrients,
              names_to = "Nutrients",
               values_to = "Nutrients_value") %>%
  mutate(Nutrients_value = ifelse(Nutrients_value <0, 0, Nutrients_value),
         Nutrients = factor(Nutrients,
                            levels = nutrients,
                            labels = c("Total\nOrganic\nNitrogen", "Nitrite", "Phosphate", "Silicate")))

ggplot(data = nuts_df,
       aes(x = as.factor(Sample_Index),
           y = Nutrients_value,
           color = Nutrients,
           fill = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(Nutrients_value, 1)), vjust = 1.5, size = 4)+
  expand_limits(y = -0.2)+
  scale_color_aaas()+
  #expand_limits(y =  7) +
 guides(fill = F)+
  xlab("Station Number")+
  ylab("Nutrients (uM)")+
  labs(title = "Nutrients - PE477 & PE486 cruises"
  )+
  facet_grid(. ~ Nutrients,
             scale = 'free')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        title = element_text(face = 'bold')
)

#Temperature and Salinity ####

ts_df <- pe_df %>%
  select(all_of(c(common, ts))) %>%
  pivot_longer(cols = ts,
               names_to = "TS",
               values_to = "TS_Value") %>%
  mutate(TS = factor(TS,
                            levels = ts,
                            labels = c("Temperature", "Salinity")))

ggplot(data = ts_df,
       aes(x = as.factor(Station_Number),
           y = TS_Value,
           color = TS,
           fill = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(TS_Value, 1)), vjust = 1.5, size = 4)+
  expand_limits(y = -2)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
 # ylab("Nutrients (uM)")+
  labs(title = "Temperature & Salinity - PE477 & PE486 cruises"
  )+
  facet_grid(Location ~ TS,
             scale = 'free')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        title = element_text(face = 'bold')
  )

#Bacteria ####
ba_df<- pe_df %>%
  select(all_of(c(common, bacteria))) %>%
  pivot_longer(cols = bacteria,
               names_to = "Bacteria",
               values_to = "Ba_abundance") %>%
  mutate(Bacteria = factor(Bacteria,
                     levels = bacteria,
                     labels = c("Total Bacteria",
                                "HNA Bacteria", "LNA Bacteria")))

ggplot(data = ba_df,
       aes(x = as.factor(Sample_Index),
           y = Ba_abundance/1e+6,
           color = as.factor(Bacteria),
           fill = Location,
           shape = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(Ba_abundance/1e+6, 1)), vjust = 2, size = 4) + 
  expand_limits(y = -0.1)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
  ylab("Bacterial abundance (in millions per mL)")+
  labs(title = "Bacterial Abundance - PE477 & PE486 cruises"
  )+
   facet_grid(. ~ Bacteria,
             scale = 'free')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
  )

#Viruses ####
vi_df<- pe_df %>%
  select(all_of(c(common, viruses))) %>%
  pivot_longer(cols = viruses,
               names_to = "Viruses",
               values_to = "Vi_abundance") %>%
  mutate(Viruses = factor(Viruses,
                           levels = viruses,
                           labels = c("Total Viruses",
                                      "V1 Viruses", "V2 Viruses",
                                      "V3 Viruses")))

ggplot(data = vi_df,
       aes(x = as.factor(Sample_Index),
           y = Vi_abundance/1e+6,
           color = as.factor(Viruses),
           fill = Location,
           shape = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(Vi_abundance/1e+6, 1)), vjust = 1.5, size = 4) + 
  expand_limits(y = -2)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
  ylab("Viral abundance (in millions per mL)")+
  labs(title = "Viral Abundance - PE477 & PE486 cruises"
  )+
  facet_grid(. ~ Viruses,
             scale = 'free')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
  )

#VBR####

ggplot(data = pe_df,
       aes(x = as.factor(Sample_Index),
           y = VBR,
           #color = as.factor(Viruses),
           color = Location,
           shape = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(VBR, 1)), vjust = 1.5, size = 4) + 
  expand_limits(y = -2)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
  ylab("VBR")+
  labs(title = "Virus bacterium ratio - PE477 & PE486 cruises"
  )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
        
  )

#VP ####

vp_df<- pe_df %>%
  select(all_of(c(common, 'VP_Lytic', 'VP_Lysogenic'))) %>%
  pivot_longer(cols = c('VP_Lytic', 'VP_Lysogenic'),
                 names_to = 'Sample_Type',
               names_prefix = 'VP_',
               values_to = 'vp_value')
 
ggplot(data = vp_df,
       aes(x = as.factor(Sample_Index),
           y = vp_value/1e+6,
           color = Location,
           shape = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(vp_value/1e+6, 3)), hjust = 1.5, size = 4, angle = 90) + 
  expand_limits(y = -0.5)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
  ylab("Viral Production Rate (in million viruses per mL per hour)")+
  labs(title = "Viral Production rate (Lytic & Lysogenic) - PE477 & PE486 cruises"
  )+
  facet_grid(. ~ Sample_Type )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
        
  )
#VP_abs ####

abs_vp_df<- pe_df %>%
  select(all_of(c(common, 'abs_VP_Lytic', 'abs_VP_Lysogenic'))) %>%
  pivot_longer(cols = c('abs_VP_Lytic', 'abs_VP_Lysogenic'),
               names_to = 'Sample_Type',
               names_prefix = 'abs_VP_',
               values_to = 'abs_vp_value')

ggplot(data = abs_vp_df,
       aes(x = as.factor(Sample_Index),
           y = abs_vp_value/1e+6,
           color = Location,
           shape = Location)) +
  geom_hline(yintercept = 0,
             size = 1,
             alpha = 0.2)+
  geom_point(size = 2)+
  geom_text(aes(label = round(abs_vp_value/1e+6, 3)), hjust = 1.5, size = 4, angle = 90) + 
  expand_limits(y = -1)+
  scale_color_aaas()+
  guides(fill = F)+
  xlab("Station Number")+
  ylab("Absolute Viral Production (in million viruses per mL)")+
  labs(title = "Absolute Viral Production (Lytic & Lysogenic) - PE477 & PE486 cruises"
  )+
  facet_grid(. ~ Sample_Type )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
        
  )

#Need to plot % of bacteria lost due to lytic and lysogenic production

#####MAPS####

#Basic plot with lats and longs first 
#Aa our statiosn are close to eachother. We'll add jitter/nudeg position for the labels.
pe_df<- pe_df %>%
  mutate(nudge_lat = Latitude - 0.1,
         nudge_long = Longitude)


pe_df[pe_df$Sample_Index == 'PE486_6',][,c('nudge_lat', 'nudge_long')]<- c(55.3, 1.55)
pe_df[pe_df$Sample_Index == 'PE486_7',][,c('nudge_lat', 'nudge_long')]<- c(55.2, 1.9)


# Get the bounding box (range) of your data
bbox <- make_bbox(lon = pe_df$Longitude, lat = pe_df$Latitude, f = 0.6)  # f is the added margin

# Retrieve the map
map_data <- get_map(location = bbox, source = "stamen", maptype = "toner-lite", zoom = 7)

# Plot the data on the map
pe_map<-ggmap(map_data) +
  geom_segment(data = pe_df  %>%
                 dplyr::filter(Sample_Index %in% c('PE486_6', 'PE486_7')), 
               aes(x = Longitude, y = Latitude, xend = nudge_long, yend = nudge_lat), 
               color = "grey30", linewidth = 0.5) +
  geom_point(data = pe_df, 
             aes(x = Longitude, y = Latitude, color = Location), size = 2.5, alpha = 0.5) +
  geom_text(data = pe_df %>%
              dplyr::filter(!Sample_Index %in% c('PE486_6', 'PE486_7')),
            aes(x = nudge_long, y = nudge_lat, label = as.character(Station_Number)), vjust = "bottom", size = 5) +
  geom_text(data = pe_df %>%
              dplyr::filter(Sample_Index %in% c('PE486_6', 'PE486_7')),
            aes(x = nudge_long, y = nudge_lat, label = as.character(Station_Number)), vjust = "bottom", size = 5) +
  scale_color_aaas()+
  labs(title = "PE477 & PE486 Sampling Stations")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
        
  )
pe_map
####Add Bathymetry from emodnet ####

bathy_raster_NorthSea <- raster::raster("Mean_depth_in_multi_colour_no_land.geotif")
bathy_agg <- aggregate(bathy_raster_NorthSea, fact = 10)
writeRaster(bathy_agg, 'NorthSea_bathymetry_raster_EMODnet.geotif', format = 'GTiff')
bathy_agg<- raster::raster("")
bathy_agg <- raster::raster("NorthSea_bathymetry_raster_EMODnet.tif")

# Get the bounding box (range) of your data
bbox <- make_bbox(lon = pe_df$Longitude, lat = pe_df$Latitude, f = 0.6)  # f is the added margin

# Retrieve the map
map_data <- get_map(location = bbox, source = "stamen", maptype = "toner-lite", zoom = 7)


ggmap(map_data) +
  geom_contour(data = as.data.frame(rasterToPoints(bathy_agg)), 
               aes(x = x,  y = y, layer = bathy_agg), color = "blue") +
  geom_segment(data = pe_df  %>%
                 dplyr::filter(Sample_Index %in% c('PE486_6', 'PE486_7')), aes(x = Longitude, y = Latitude, xend = nudge_long, yend = nudge_lat), color = "grey30", size = 0.5) +
  geom_point(data = pe_df, 
             aes(x = Longitude, y = Latitude, color = Location), size = 2.5, alpha = 0.5) +
  geom_text(data = pe_df %>%
              dplyr::filter(!Sample_Index %in% c('PE486_6', 'PE486_7')),
            aes(x = nudge_long, y = nudge_lat, label = as.character(Station_Number)), vjust = "bottom", size = 3) +
  geom_text(data = pe_df %>%
              dplyr::filter(Sample_Index %in% c('PE486_6', 'PE486_7')),
            aes(x = nudge_long, y = nudge_lat, label = as.character(Station_Number)), vjust = "bottom", size = 3) +
  scale_color_aaas()+
  labs(title = "PE477 & PE486 Sampling Stations")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = "#202a47"),
        strip.text = element_text(color = 'white',
                                  size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')
        
  )

#cheecking if they have the saem coordinate systems
crs(bathy_raster_NorthSea)
st_crs()
#Adding grobs for different variables on map

#we'll have to adjust thenudge position accordingly later.

#for nutrients, a barchart with values on them would be good
#lets subset one station to get strted
sdf_nuts<- nuts_df %>%
  dplyr:: filter(Sample_Index == 'PE477_1') 


ggplot(data = sdf_nuts,
       aes(x = Nutrients,
           y = Nutrients_value,
           fill = Nutrients)) +
  geom_bar(stat = 'identity')+
  scale_fill_aaas()+
  ylab("value (uM)")+
  theme_classic()+
  theme(legend.position = 'none')

unique_samples <- unique(nuts_df$Sample_Index)
plots_nuts <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(nuts_df$Nutrients_value, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  nut <- unique_samples[i]
  
  nuts_plot <- ggplot(data = nuts_df %>% 
                        dplyr::filter(Sample_Index == nut),
                      aes(x = Nutrients,
                          y = Nutrients_value,
                          fill = Nutrients)) +
    geom_bar(stat = 'identity',
             width = 0.7,
             position = position_dodge(width = 0.05)) +
    geom_text(aes(label = sprintf("%.2f", Nutrients_value)), vjust = -0.5, size = 3) +  # Display values with 2 decimal places
    scale_fill_aaas() +
   # ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_classic() +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank())+
    coord_fixed(ratio = 0.5)
  
  plots_nuts[[i]] <- nuts_plot
}
plots_nuts[[7]]

# Assuming you've created bbox using the make_bbox() function
bbox <- make_bbox(lon = pe_df$Longitude, lat = pe_df$Latitude, f = 0.6)

# Convert the numeric vector to a named list
bbox <- list(
  left = bbox[1],
  bottom = bbox[2],
  right = bbox[3],
  top = bbox[4]
)

# Function to add an inset plot to a given map
add_plot_to_map <- function(main_plot, inset_plot, x, y, width = 0.2, height = 0.2) {
  print(paste("Placing plot at coordinates:", x, y))
  vp <- grid::viewport(width = width, height = height, x = x, y = y, just = c("center", "center"))
  print(inset_plot, vp = vp)
  return(main_plot)
}


# Updated script to overlay each nutrient plot
plot_data <- nuts_df %>%
  distinct(Sample_Index) %>%
  left_join(pe_df, by = "Sample_Index") %>%
  mutate(Plot = plots_nuts)

# Overlay each nutrient plot on the main map
final_nuts_plot <- purrr::reduce(
  seq_along(plots_nuts),
  function(p, i) {
    add_plot_to_map(p, plot_data$Plot[[i]], plot_data$nudge_long[i], plot_data$nudge_lat[i], width = 1, height = 1)
  },
  .init = pe_map
)

print(final_nuts_plot)


#PCA####

subset_pe_df <- pe_df %>%
  select(-all_of(c('Location', 'Station_Number',
  'Depth', 'Time_Range', 'Population',
                   'Bacterial_Sample_Name', 'Viral_Sample_Name',
                   'Expt_Date', 'Latitude', 'Longitude',
                   'nudge_lat', 'nudge_long', 'Sample_Index',
  'Temperature', 'Salinity', 'Nitrate', 'Total_Bacteria', 'Total_Viruses', 'VBR'))) 

cols_to_check <- c("TON", "Nitrite", "Phosphate", "Silicate")

# Check for negative values and replace with 0
for (col in cols_to_check) {
  subset_pe_df[[col]][subset_pe_df[[col]] < 0] <- 0
}

#to handle missing values
# Impute the data
imputed_data <- mice(subset_pe_df, m=5, maxit=50, method='pmm', seed=500)
completed_data <- complete(imputed_data, action = 1)
scaled_data <- scale(subset_pe_df)
pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

# Scree Plot
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
barplot(explained_var, main = "Scree Plot", 
        xlab = "Principal Component", ylab = "Percentage of Variance Explained",
        ylim = c(0, max(explained_var) + 10), col = "steelblue")
biplot(pca_result, main = "PCA Biplot", cex = 0.7)

scores <- data.frame(pca_result$x[, 1:2])
scores$Sample_Index <- pe_df$Sample_Index
scores$Location<- pe_df$Location
scores$Station_Number<- pe_df$Station_Number

loadings <- data.frame(pca_result$rotation[, 1:2])


ggplot(scores, aes(x = PC1, y = PC2, color = Location)) +
  geom_segment(aes(x = 0, y = 0, xend = 10*PC1, yend = 10*PC2), data = loadings, arrow = arrow(type = "closed", length = unit(0.2, "cm")), alpha = 0.5, color = "black") +
  geom_text(data = loadings, aes(x = 10*PC1, y = 10*PC2, label = rownames(loadings)), color = "black")+
  geom_point(aes(), size = 3) + 
  scale_color_aaas()+
  geom_text(aes(label = Station_Number), vjust = -0.5, hjust = 1) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_bw()


#Downloading subplots to plot on INkscape ####

#Nutrients - together ####

sdf_nuts<- nuts_df %>%
  dplyr:: filter(Sample_Index == 'PE477-1') 


ggplot(data = sdf_nuts,
       aes(x = Nutrients,
           y = Nutrients_value,
           fill = Nutrients)) +
  geom_bar(stat = 'identity')+
  scale_fill_aaas()+
  ylab("value (uM)")+
  theme_classic()+
  theme(legend.position = 'right')

unique_samples <- unique(nuts_df$Sample_Index)
unique_nutrients<- unique(nuts_df$Nutrients)
plots_nuts <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(nuts_df$Nutrients_value, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  for (j in seq_along(unique_nutrients)) {
    nut <- unique_samples[i]
    nutrient <- unique_nutrients[j]
  
  nuts_plot <- ggplot(data = nuts_df %>% 
                        dplyr::filter(Sample_Index == nut &
                                        Nutrients == nutrient),
                      aes(x = Nutrients,
                          y = Nutrients_value,
                          fill = Nutrients)) +
    geom_bar(stat = 'identity',
             width = 0.7,
             position = position_dodge(width = 0.05)) +
    geom_text(aes(label = sprintf("%.2f", Nutrients_value)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    # ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'none',
          plot.caption.position = "plot")+
    coord_fixed(ratio = 0.5)
  
  plots_nuts[[i]] <- nuts_plot
  
  filename <- paste0("figures/nuts_subplots/", nut, ".svg")  # Construct the filename using the current nut value
  ggsave(filename, plot = nuts_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  }
}

#Nutrients separte ####


#Nutrients - together ####
nuts_df <- pe_df %>%
  select(all_of(c(common, nutrients))) %>%
  pivot_longer(cols = nutrients,
               names_to = "Nutrients",
               values_to = "Nutrients_value") %>%
  mutate(Nutrients_value = ifelse(Nutrients_value <0, 0, Nutrients_value),
         Nutrients = factor(Nutrients,
                            levels = nutrients,
                            labels = c("TON", "Nitrite", "Phosphate", "Silicate")))

unique_samples <- unique(nuts_df$Sample_Index)
unique_nutrients<- unique(nuts_df$Nutrients)
plots_nuts <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(nuts_df$Nutrients_value, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  for (j in seq_along(unique_nutrients)) {
    nut <- unique_samples[i]
    nutrient <- unique_nutrients[j]
    
    nuts_plot <- ggplot(data = nuts_df %>% 
                          dplyr::filter(Sample_Index == nut &
                                          Nutrients == nutrient) 
                        ,
                        aes(x = Nutrients,
                            y = Nutrients_value,
                            fill = Nutrients)) +
      geom_bar(stat = 'identity',
               width = 1.5,
               position = position_dodge(width = 0.05)) +
      geom_text(aes(label = sprintf("%.2f", Nutrients_value)), 
                vjust = -0.5, size = 3,
                fontface = 'bold') +  # Display values with 2 decimal places
     labs(caption = unique_samples[i])+
      scale_fill_aaas() +
      # ylab("value (uM)") +
      ylim(0, global_ymax + 1) + 
      theme_void() +
      theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
            legend.position = 'right',
            plot.caption.position = "plot")+
      coord_fixed(ratio = 0.5)
    
    plots_nuts[[i]] <- nuts_plot
    
    filename <- paste0("figures/nuts_subplots/", nut, "_", nutrient, ".svg")  # Construct the filename using the current nut value
    #ggsave(filename, plot = nuts_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  }
}


#Temperature and salinity. individual graphs####

ts_df
unique_samples <- unique(ts_df$Sample_Index)
unique_TS<- unique(ts_df$TS)
plots_ts <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(ts_df$TS_Value, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  for (j in seq_along(unique_TS)) {
    sample <- unique_samples[i]
    ts <- unique_TS[j]
    
    ts_plot <- ggplot(data = ts_df %>% 
                          dplyr::filter(Sample_Index == sample &
                                          TS == ts) 
                        ,
                        aes(x = TS,
                            y = TS_Value,
                            fill = TS)) +
      geom_bar(stat = 'identity',
               width = 1,
               position = position_dodge(width = 0.05)) +
      geom_text(aes(label = sprintf("%.2f", TS_Value)), 
                vjust = -0.5, size = 3,
                fontface = 'bold') +  # Display values with 2 decimal places
      # labs(caption = unique_samples[i])+
      scale_fill_aaas() +
      ylab("value (uM)") +
      ylim(0, global_ymax + 1) + 
      theme_void() +
      theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
            legend.position = 'none',
            plot.caption.position = "plot")+
      coord_fixed(ratio = 0.5)
    
    plots_ts[[i]] <- ts_plot
    
    filename <- paste0("figures/ts_subplots/", sample, "_", ts, ".svg")  # Construct the filename using the current nut value
    ggsave(filename, plot = ts_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  }
} #DOING THIS LATER

ts_df
unique_samples <- unique(ts_df$Sample_Index)
unique_TS<- unique(ts_df$TS)
plots_ts <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(ts_df$TS_Value, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
    sample <- unique_samples[i]
    
    
    ts_plot <- ggplot(data = ts_df %>% 
                        dplyr::filter(Sample_Index == sample),
                      aes(x = TS,
                          y = TS_Value,
                          fill = TS)) +
      geom_bar(stat = 'identity',
               width = 0.7,
               position = position_dodge(width = 0.1)
               ) +
      geom_text(aes(label = sprintf("%.2f", TS_Value)), 
                vjust = -0.5, size = 3,
                fontface = 'bold') +  # Display values with 2 decimal places
      labs(caption = unique_samples[i])+
      scale_fill_aaas() +
      ylab("value (uM)") +
      ylim(0, global_ymax + 1) + 
      theme_void() +
      theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
            legend.position = 'none',
            plot.caption.position = "plot")
    
    plots_ts[[i]] <- ts_plot
    
    filename <- paste0("figures/ts_subplots/", sample, ".svg")  # Construct the filename using the current nut value
    ggsave(filename, plot = ts_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 

}

#Bacteria #### map

ba_df
unique_samples <- unique(ba_df$Sample_Index)
unique_Bacteria<- unique(ba_df$Bacteria)
plots_ba <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(ba_df$Ba_abundance/1e+6, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
  sample <- unique_samples[i]
  
  
  ba_plot <- ggplot(data = ba_df %>% 
                      dplyr::filter(Sample_Index == sample),
                    aes(x = Bacteria,
                        y = Ba_abundance/1e+6,
                        fill = Bacteria)) +
    geom_bar(stat = 'identity',
             width = 0.7,
             position = position_dodge(width = 0.1)
    ) +
    geom_text(aes(label = sprintf("%.2f", Ba_abundance/1e+6)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'none',
          plot.caption.position = "plot")
  
  plots_ba[[i]] <- ba_plot
  
  filename <- paste0("figures/ba_subplots/", sample, ".svg")  # Construct the filename using the current nut value
  ggsave(filename, plot = ba_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  
}

#Viruses #### map

vi_df
unique_samples <- unique(vi_df$Sample_Index)
unique_Viruses<- unique(vi_df$Viruses)
plots_vi <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(vi_df$Vi_abundance/1e+6, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
  sample <- unique_samples[i]
  
  
  vi_plot <- ggplot(data = vi_df %>% 
                      dplyr::filter(Sample_Index == sample),
                    aes(x = Viruses,
                        y = Vi_abundance/1e+6,
                        fill = Viruses)) +
    geom_bar(stat = 'identity',
             width = 0.5,
             position = position_dodge(width = 0.0)
    ) +
    geom_text(aes(label = sprintf("%.2f", Vi_abundance/1e+6)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'none',
          plot.caption.position = "plot")
  
  plots_vi[[i]] <- vi_plot
  
  filename <- paste0("figures/vi_subplots/", sample, ".svg")  # Construct the filename using the current nut value
  ggsave(filename, plot = vi_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  
}


#VBR####

pe_df
unique_samples <- unique(pe_df$Sample_Index)
plots_vbr <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(pe_df$VBR, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
  sample <- unique_samples[i]
  
  
  vbr_plot <- ggplot(data = pe_df %>% 
                      dplyr::filter(Sample_Index == sample),
                    aes(x = 1,
                        y = VBR)) +
    geom_bar(stat = 'identity',
             aes(width = 0.1),
             position = position_dodge(width = 0.0),
             fill = 'darkgreen'
    ) +
    geom_text(aes(label = sprintf("%.2f", VBR)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'none',
          plot.caption.position = "plot")
  
  plots_vbr[[i]] <- vbr_plot
  
  filename <- paste0("figures/vbr_subplots/", sample, ".svg")  # Construct the filename using the current nut value
  ggsave(filename, plot = vbr_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  
}

## VP #### map
vp_df<- pe_df %>%
  select(all_of(c(common, 'VP_Lytic', 'VP_Lysogenic'))) %>%
  pivot_longer(cols = c('VP_Lytic', 'VP_Lysogenic'),
               names_to = 'Sample_Type',
               names_prefix = 'VP_',
               values_to = 'vp_value')

unique_samples <- unique(vp_df$Sample_Index)

plots_vp <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(vp_df$vp_value/1e+3, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
  sample <- unique_samples[i]
  
  
  vp_plot <- ggplot(data = vp_df %>% 
                      dplyr::filter(Sample_Index == sample),
                    aes(x = Sample_Type,
                        y = vp_value/1e+3,
                        fill = Sample_Type)) +
    geom_bar(stat = 'identity',
             width = 0.5,
             position = position_dodge(width = 0.0)
    ) +
    geom_text(aes(label = sprintf("%.2f", vp_value/1e+3)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    ylab("value (uM)") +
    ylim(0, global_ymax + 1) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'right',
          plot.caption.position = "plot")
  
  plots_vp[[i]] <- vp_plot
  
  filename <- paste0("figures/vp_subplots/", sample, ".svg")  # Construct the filename using the current nut value
  #ggsave(filename, plot = vp_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  
}


## ABS VP #### map
abs_vp_df<- pe_df %>%
  select(all_of(c(common, 'VP_Lytic', 'VP_Lysogenic'))) %>%
  pivot_longer(cols = c('VP_Lytic', 'VP_Lysogenic'),
               names_to = 'Sample_Type',
               names_prefix = 'VP_',
               values_to = 'abs_vp_value')

unique_samples <- unique(abs_vp_df$Sample_Index)

plots_abs_vp <- vector("list", length(unique_samples))

# Determine the global y-axis maximum limit
global_ymax <- max(abs_vp_df$abs_vp_value/1e+3, na.rm = TRUE)

for (i in seq_along(unique_samples)) {
  
  sample <- unique_samples[i]
  
  
  abs_vp_plot <- ggplot(data = abs_vp_df %>% 
                      dplyr::filter(Sample_Index == sample),
                    aes(x = Sample_Type,
                        y = abs_vp_value/1e+3,
                        fill = Sample_Type)) +
    geom_bar(stat = 'identity',
             width = 0.5,
             position = position_dodge(width = 0.0)
    ) +
    geom_text(aes(label = sprintf("%.2f", abs_vp_value/1e+3)), 
              vjust = -0.5, size = 3,
              fontface = 'bold') +  # Display values with 2 decimal places
    labs(caption = unique_samples[i])+
    scale_fill_aaas() +
    ylab("value (uM)") +
    ylim(0, global_ymax + 10) + 
    theme_void() +
    theme(plot.caption = element_text(hjust = 0.5, face = "bold"),  
          legend.position = 'none',
          plot.caption.position = "plot")
  
  plots_abs_vp[[i]] <- abs_vp_plot
  
  filename <- paste0("figures/abs_vp_subplots/", sample, ".svg")  # Construct the filename using the current nut value
  ggsave(filename, plot = abs_vp_plot, width = 5, height = 5, dpi = 96, unit = 'cm') 
  
}

#abs_VP vs VP

global_min <- min(min(pe_df$VP_Lytic, na.rm = TRUE), min(pe_df$abs_VP_Lytic, na.rm = TRUE))
global_max <- max(max(pe_df$VP_Lytic, na.rm = TRUE), max(pe_df$abs_VP_Lytic, na.rm = TRUE))

vp_rate<- ggplot(data  = pe_df,
       aes(x = VP_Lytic*24,
       y = abs_VP_Lytic)) +
  geom_point()+
  theme_bw()+
  geom_abline(slope = 1)+
  geom_smooth(method = 'lm', color = 'maroon')+
  xlim(global_min, global_max) +
  ylim(global_min, global_max)

vp_abs<- ggplot(data  = pe_df,
       aes(x = VP_Lysogenic*24,
           y = abs_VP_Lysogenic)) +
  geom_point()+
  theme_bw()+
  geom_abline(slope = 1)+
  geom_smooth(method = 'lm', color = 'maroon')+
  xlim(global_min, global_max) +
  ylim(global_min, global_max)

vp_rate + vp_abs


#Percent of BACTERIA LYSED



  
