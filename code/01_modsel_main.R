library(tidyverse)
library(terra)
library(tidyterra)
library(car)
library(lme4)
library(parallel)

set.seed(64392)

source("code/modsel_helper.R")

#### Setup ####

covars <- c("Longitude", "Latitude", "NDVI", "NDWI", "Elevation", "Slope",
            "Percent_Deciduous", "Percent_Evergreen", "Percent_Mixed",
            "Percent_Impervious", "Percent_Wetlands")
scales <- c(25, 250)

datlist <- readRDS("data/bobcat_datlist.RDS")

glmm_result <- do_modsel(datlist, ranef = "(1 | siteID)", scales, ncores = 10)
saveRDS(glmm_result, "glmm_result.RDS")



