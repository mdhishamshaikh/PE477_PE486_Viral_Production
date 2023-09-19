library(tidyverse)
library(leaflet)
library(ggsci)
library(ggmap)
library(sf)
library(cowplot)

pe_df<- read.csv("./PE_Cruises/results/PE_Cruises_viral_production/vp_abundance_nutrients_ts.csv")
pe_df <- pe_df %>%
  mutate(Sample_Index = factor(paste0(Location, "_", Station_Number)))

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
  mutate(Bacteria = factor(Viruses,
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

