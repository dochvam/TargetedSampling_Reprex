library(unmarked)
library(lme4)
library(tidyverse)
library(terra)
library(tidyterra)
library(gridExtra)

# mapview::mapview(all_pts)

spatial_buffer <- 250
n_targeted_cameras <- 20
targeted_pts_distance <- 75

source("code/modsel_helper.R")

glmm_reslist <- readRDS("glmm_result.RDS")

datlist <- readRDS("data/bobcat_datlist.RDS")

best_model_df <- glmm_reslist$best_mod_df

covars_base <- c("NDVI", "NDWI", "Elevation", "Slope",
                 "Percent_Deciduous", "Percent_Evergreen", "Percent_Mixed",
                 "Percent_Impervious")
combos <- expand.grid(covars_base, c("_25", "_250"))
covars <- c("Longitude", "Latitude", paste0(combos[, 1], combos[, 2]))
covars_base  <- c("Longitude", "Latitude", covars_base)

formula <- formula_from_df_glmm(best_model_df, "(1|siteID)")

best_fit <- glmer(formula, datlist[[1]], family = "binomial", control = glmerControl(optimizer = "bobyqa"))

summary(best_fit)
car::vif(best_fit)



#### Make plots of covariate predictions ####
newdata_baseline <- datlist[[1]][1:20, c(covars, paste0(covars, "_sq"))]
newdata_baseline[] <- 0

panels <- list()
best_model_df <- best_model_df %>% filter(state > 0)
for (i in 1:nrow(best_model_df)) {
  this_cov <- best_model_df$param[i]
  if (!is.na(best_model_df$scale[i])) this_cov <- paste0(this_cov, "_", best_model_df$scale[i])
  
  this_sd <- datlist$scaling_factors$sd[datlist$scaling_factors$covar == this_cov]
  this_mean <- datlist$scaling_factors$mean[datlist$scaling_factors$covar == this_cov]
  submod <- best_model_df$submodel[i]
  
  this_range <- range(datlist[[1]][, this_cov])
  
  newdata <- newdata_baseline
  newdata[, this_cov] <- seq(from = this_range[1], to = this_range[2], length.out = 20)
  newdata[, paste0(this_cov, "_sq")] <- seq(from = this_range[1], to = this_range[2], length.out = 20)^2
  
  
  predicted <- predict(best_fit,
                       newdata = newdata, 
                       type = "link", 
                       re.form = NA, se.fit = TRUE) %>% 
    as.data.frame() %>% 
    mutate(covar_val = unlist(newdata[, this_cov])) %>% 
    mutate(covar_val_unscaled = covar_val * this_sd + this_mean) %>%
    mutate(
      Predicted = plogis(fit),
      lower = plogis(fit - 1.96*se.fit), 
      upper = plogis(fit + 1.96*se.fit)
    )
  
  panels[[i]] <- ggplot(predicted, aes(covar_val_unscaled, Predicted, ymin = lower, ymax = upper)) +
    geom_line() +
    geom_ribbon(alpha = 0.2) +
    theme_minimal() +
    xlab(this_cov) + ylab("Observation prob.")
}

covar_plot <- arrangeGrob(grobs = panels, nrow = 4)
ggsave("plots/covar_effects_plot.jpg", covar_plot, width = 8, height = 10)

#### Predict across the whole study region ####

covar_brick_sm <- rast("data/covar_brick_AOI.tif")

south_mountains <- vect("data/SouthMtnsParkBoundary.shp") %>% 
  project(crs(covar_brick_sm))


newdata <- covar_brick_sm %>% as.matrix() %>% as.data.frame() %>% na.omit()

longlat <- as.data.frame(terra::project(vect(crds(covar_brick_sm), crs = crs(covar_brick_sm)), "EPSG:4326"), geom = "XY")
colnames(longlat) <- c("Longitude", "Latitude")

newdata <- newdata %>% 
  bind_cols(longlat)

scaling_factors <- datlist$scaling_factors
for (i in 1:nrow(scaling_factors)) {
  newdata[, scaling_factors$covar[i]] <- 
    (newdata[, scaling_factors$covar[i]] - scaling_factors$mean[i]) /
    scaling_factors$sd[i]
  
  newdata[, paste0(scaling_factors$covar[i], "_sq")] <- newdata[, scaling_factors$covar[i]]^2
}

predicted <- predict(best_fit, type = "response", newdata = newdata, re.form = NA)

predicted_link <- predict(best_fit, type = "link", newdata = newdata, re.form = NA, se.fit = TRUE)

predicted_link$fit
predicted_link$se

newdata$point_pred <-  predicted


predicted_raster <- covar_brick_sm[[1]]
predicted_raster_link <- covar_brick_sm[[1]]
predicted_raster_se <- covar_brick_sm[[1]]

values(predicted_raster)[!is.na(values(predicted_raster))] <- newdata$point_pred

values(predicted_raster_link)[!is.na(values(predicted_raster_link))] <- unlist(predicted_link$fit)
values(predicted_raster_se)[!is.na(values(predicted_raster_se))] <- unlist(predicted_link$se)

snapshot_pts <- read_csv("data/SM_2024_full_grid.csv") %>% 
  mutate(type = "Systematic") %>% 
  rename(Longitude = Long, Latitude = Lat) %>% 
  vect(geom = c("Longitude", "Latitude"),
       crs = "+proj=longlat") %>% 
  project(crs(predicted_raster))

predicted_raster_of <- predicted_raster_link
values(predicted_raster_of) <- plogis(values(predicted_raster_link) - 1.96*values(predicted_raster_se))

pred_plot <- ggplot() +
  geom_spatraster(data = predicted_raster) +
  geom_spatvector(data = snapshot_pts, col = "orangered", pch = 18, size = 3) +
  scale_fill_viridis_c() +
  ggtitle("Predicted prob. of detecting a bobcat on a 10-day survey")
ggsave("plots/detrate_pred_plot.jpg", pred_plot, width = 5, height = 4)


plot_all_covars <- FALSE # Don't need to do this every time
if (plot_all_covars) {
  covar_plots <- list()
  for (i in 1:nlyr(covar_brick_sm)) {
    this_name <- names(covar_brick_sm)[i]
    covar_plots[[i]] <- ggplot() + 
      geom_spatraster(data = covar_brick_sm[[i]]) +
      scale_fill_viridis_c(this_name)
    ggsave(paste0("plots/covar_plots/", this_name, ".jpg"),
           width = 6, height = 4, dpi = 200)
  }
}


#### Choose the targeted points! ####
targeted_points_df <- data.frame(x = numeric(0), y = numeric(0))
temp_pred_ras <- predicted_raster_of
for (i in 1:n_targeted_cameras) {
  # Find the point with the highest predicted value
  best_cell <- which(values(temp_pred_ras) == max(values(temp_pred_ras), na.rm = T))
  best_cell <- best_cell[1]
  
  coord <- xyFromCell(temp_pred_ras, best_cell)
  buf <- vect(coord) %>% buffer(width = 75)
  cells_in_buffer <- cells(temp_pred_ras, buf)
  cell_ids_within_d <- cells_in_buffer[,2]
  
  # adj_cells <- adjacent(temp_pred_ras, best_cell, directions = 8)
  
  values(temp_pred_ras)[best_cell] <- 0
  values(temp_pred_ras)[cell_ids_within_d] <- 0
  
  targeted_points_df <- rbind(targeted_points_df,
                              coord)
}

targeted_pts <- targeted_points_df %>% 
  mutate(type = "Targeted") %>% 
  vect(geom = c("x", "y"), crs = crs(predicted_raster))

if (nrow(snapshot_pts) < n_targeted_cameras) {
  inds <- rep(1:nrow(snapshot_pts), 10)
  snapshot_pts <- snapshot_pts[inds[1:n_targeted_cameras]]
}

all_pts <-  rbind(targeted_pts, snapshot_pts)

#### Plot systematic vs. targeted predictions ####

panel_a <- ggplot() +
  geom_spatraster(data = predicted_raster) +
  geom_spatvector(data = all_pts, 
                  aes(color = type), size = 3, pch = 18) +
  scale_fill_viridis_c("Obs. prob.", option = "mako") +
  # ggtitle("Predicted prob. of detecting a bobcat on a 10-day survey") +
  scale_color_manual("Array type", values = c("#E09F3E", "#D33F49")) +
  theme_minimal() +
  ggtitle("A. Systematic and targeted array placements")
ggsave("plots/detrate_pred_plot.jpg", pred_plot, width = 5, height = 4)



all_pts$pred_rate <- extract(predicted_raster, all_pts)[, 2]
all_pts$al_one <- 1 - dbinom(0, 6, all_pts$pred_rate)

writeVector(all_pts, "data/SM_sampling_locs_GLMM.shp", overwrite = TRUE)


panel_b <- ggplot(as.data.frame(all_pts)) + 
  geom_histogram(aes(pred_rate, fill = type, group = type)) +
  scale_fill_manual("Array type", values = c("#E09F3E", "#D33F49")) +
  theme_minimal() +
  xlab("Prob. of at least one detection per 10-day window") + ylab("Number of cameras") +
  ggtitle("B. Distribution of predicted det. probs.")

all_pts %>% 
  as.data.frame() %>% 
  group_by(type) %>% 
  summarize(mean = mean(al_one))
# 
# ggplot() +
#   geom_spatvector(data = all_pts[all_pts$type == "Systematic",], aes(color = predicted_rate), size = 5) +
#   geom_spatvector_label(data = all_pts[all_pts$type == "Systematic",], aes(label = round(predicted_rate, 2)), size = 3)
#   

# mapview(all_pts, zcol = "type")


# Expected distribution of number of 7-day intervals containing bobcat,
# if every camera runs for 42 days

systematic_array_probs <- rep(all_pts$pred_rate[all_pts$type == "Systematic"], 6)
targeted_array_probs   <- rep(all_pts$pred_rate[all_pts$type == "Targeted"], 6)

pred_df <- data.frame(
  count = 0:(14*5),
  Targeted = NA,
  Systematic = NA
)


for (i in 1:nrow(pred_df)) {
  pred_df$Targeted[i] <- poisbinom::dpoisbinom(pred_df$count[i], pp = targeted_array_probs)
  pred_df$Systematic[i] <- poisbinom::dpoisbinom(pred_df$count[i], pp = systematic_array_probs)
}

panel_c <- pred_df %>% 
  pivot_longer(cols = c(Targeted, Systematic)) %>% 
  ggplot() +
  geom_line(aes(count, value, col = name, group = name)) +
  scale_color_manual("Array type", values = c("#E09F3E", "#D33F49")) +
  theme_minimal() +
  xlim(0, 40) +
  xlab("Total 7-day intervals with detections") +
  ylab("Probability") + 
  ggtitle("C. Predicted distribution of total observations")


layout_mtx <- matrix(c(1,1,2,3),nrow=2)
plot_all <- arrangeGrob(panel_a, panel_b, panel_c, layout_matrix = layout_mtx, widths = c(1.2, 1))

ggsave(plot = plot_all, filename = "plots/sampling_plan.jpg", width = 12, height = 6, dpi = 200)
