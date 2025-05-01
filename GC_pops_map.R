# Map for Growth Chamber data
library(dplyr)
library(geodata)
library(sf)
library(terra)
library(ggspatial)
library(ggplot2)
library(ggpubr)
library(grid)
library(usmap)

setwd("C:/Users/derek/OneDrive/Documents/Boechera/Growth_Chambers/SICB_Manuscript/")
pops<- read.table("GC_Pops.csv", header = T)
pops_sf <- st_as_sf(pops, coords = c("Longitude", "Latitude"), crs = 4326)

geodata_path(getwd()) # Put elevation data in working directory
e <- ext(-107.2, -106.7, 38.7, 39.1)
srtm_data <- elevation_3s(lon=-107, lat=39)
srtm_crop <- crop(srtm_data, e)

plot(srtm_crop)

# Convert raster to a data frame for ggplot
srtm_df <- as.data.frame(srtm_crop, xy = TRUE)
colnames(srtm_df) <- c("long", "lat", "elevation")

custom_colors <- c("white", "gray30")

# New map to adjust for gridarrange
rmbl <- ggplot() +
  geom_raster(data = srtm_df, aes(x = long, y = lat, fill = elevation)) +
  scale_fill_gradientn(colors = rev(custom_colors)) +
  #scale_fill_viridis_c()+
  geom_sf(data = pops_sf,  size =1, stroke =1) +
  labs(fill = "Elevation (m)") +
  annotation_scale(location="bl", pad_x = unit(0.6, "cm"), pad_y = unit(0.63, "cm"), style = "bar")+
  theme_void()+
  coord_sf() +
  theme(legend.position = "bottom",  # Put elevation legend at the bottom
        axis.title = element_blank(),
        axis.text = element_text(size = 12)) 
rmbl

# US Map
extent_coords <- data.frame( lon = c(-107.2, -106.7), lat = c(38.7, 39.2))
extent_map <- usmap_transform(extent_coords)

# Gothic as midpoint for CO point on US map
g <- pops_sf %>% filter(Population == "283") 

# Nest CO points in US map
US_with_CO <- plot_usmap(exclude=c("AK", "HI"), color ="gray80", fill = "gray98") +
  geom_sf(data = st_geometry(g), shape = 18, size = 3, color = "gray50", inherit.aes = FALSE)
US_with_CO


# Nest US map within RMBL in the top-right corner
rmbl_with_us <- rmbl +
  annotation_custom(
    ggplotGrob(US_with_CO),
    xmin = -106.9, xmax = -106.7, 
    ymin = 39, ymax = 39.1   
  )

# Arrange plots in a single row
figure <- ggarrange(
  rmbl_with_us,
  nrow = 1,      
  widths = c(3, 1) 
)


figure
