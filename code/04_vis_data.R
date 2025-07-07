library(terra)
library(tidyverse)
library(tidyterra)
library(ggspatial)
library(cowplot)


covar_brick_sm <- rast("data/covar_brick_AOI.tif")

park_boundary <- vect("data/SouthMtnsParkBoundary.shp") %>% 
  project("EPSG:3857")

nc_map <- geodata::gadm("USA", level = 1, resolution = 1, path = "data") %>% 
  filter(NAME_1 == "North Carolina") %>% 
  project("EPSG:3857")

systematic_points <- read_csv("data/SM_2024_full_grid.csv") %>% 
  vect(geom = c("Long", "Lat"), crs = "+proj=longlat") %>% 
  project(crs(covar_brick_sm))

rast_temp <- covar_brick_sm[[2]]
values(rast_temp) <- 0
study_rect <- rast_temp %>% 
  as.polygons()
crs(study_rect) <- crs(covar_brick_sm)
study_rect <- study_rect %>% 
  project("EPSG:3857")


sm_centroid <- crds(centroids(park_boundary))
label1_df <- data.frame(
  x = sm_centroid[1],
  y = sm_centroid[2] + 22000,
  label = "South Mountains SP"
)

nc_plot <- ggplot() +
  annotation_map_tile(type = "hotstyle", zoomin = 1) +  # Light terrain
  geom_spatvector(data = nc_map, fill = "NA", color = "#333333", linewidth = 0.8) +
  geom_spatvector(data = park_boundary, color = "black", fill = "black") +
  geom_spatial_text(data = label1_df, aes(x = x, y = y, label = label),
            size = 3.3, fontface = "bold", color = "black", crs = crs(park_boundary)) +
  theme_minimal() +
  annotation_scale(
    location = "bl",     # "bl" = bottom left, options: "tl", "tr", "br"
    width_hint = 0.3     # How wide the scale bar is (relative to plot width)
  )+ xlab("") + ylab("")


label2_df <- data.frame(
  x = sm_centroid[1],
  y = sm_centroid[2],
  label = "South Mountains SP"
)
sr_centroid <- crds(centroids(study_rect))
label3_df <- data.frame(
  x = sr_centroid[1] - 500,
  y = sr_centroid[2] + 1410,
  label = "Study region"
)

sm_plot <- ggplot() +
  annotation_map_tile(type = "hotstyle", zoomin = 0) +  # Light terrain
  geom_spatvector(data = park_boundary, color = "black", fill = NA, linewidth = 1) +
  geom_spatvector(data = study_rect, fill = NA, color = "darkblue", linewidth = 1) +
  geom_spatial_text(data = label2_df, aes(x = x, y = y, label = label),
                    size = 4.5, fontface = "bold", color = "black", crs = crs(park_boundary)) +
  geom_spatial_text(data = label3_df, aes(x = x, y = y, label = label),
                    size = 2.7, fontface = "bold", color = "darkblue", crs = crs(park_boundary)) +
  annotation_scale(
    location = "bl",     # "bl" = bottom left, options: "tl", "tr", "br"
    width_hint = 0.3     # How wide the scale bar is (relative to plot width)
  ) + xlab("") + ylab("")


map_arr <- cowplot::plot_grid(nc_plot, sm_plot, ncol =1, align = "v", rel_heights = c(0.66, 1))
ggsave("plots/study_area.jpg", map_arr, width = 6, height = 7, dpi = 600)
