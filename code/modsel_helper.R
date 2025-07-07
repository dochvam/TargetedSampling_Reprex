library(unmarked)
library(tidyverse)
library(terra)
library(tidyterra)
library(car)
library(lme4)
library(parallel)

vif_threshold <- 5


formula_from_df_glmm <- function(df, ranef) {
  df$parname <- NA
  for (i in 1:nrow(df)) {
    if (df$state[i] == 1) {
      if (!is.na(df$scale[i])) {
        df$parname[i] <- paste0(df$param[i], "_", df$scale[i])
      } else {
        df$parname[i] <- paste0(df$param[i])
      }
    } else if (df$state[i] == 2) {
      if (!is.na(df$scale[i])) {
        df$parname[i] <- paste0(df$param[i], "_", df$scale[i], " + ", 
                                df$param[i], "_", df$scale[i], "_sq")
      } else {
        df$parname[i] <- paste0(df$param[i], " + ", 
                                df$param[i], "_sq")
      }
    }
  }
  
  if (sum(!is.na(df$parname)) == 0) {
    this_formula <- paste0("y ~ ", ranef)
  } else {
    this_formula <- paste0("y ~ ", ranef, " + ",
                           paste(df$parname[!is.na(df$parname)], collapse = " + "))
  }
  
  as.formula(this_formula)
}

make_proposals_glmm <- function(df, levels = c(25, 250)) {
  proposal_list <- list()
  ct <- 0
  
  for (i in 1:nrow(df)) {
    if (df$state[i] == 0) {
      if (!df$param[i] %in% c("Longitude", "Latitude")) {
        for (j in 1:length(levels)) {
          ct <- ct + 1
          proposal_list[[ct]] <- df
          proposal_list[[ct]]$state[i] <- 1
          proposal_list[[ct]]$scale[i] <- levels[j]
          ct <- ct + 1
          proposal_list[[ct]] <- df
          proposal_list[[ct]]$state[i] <- 2
          proposal_list[[ct]]$scale[i] <- levels[j]
        }
      } else {
        ct <- ct + 1
        proposal_list[[ct]] <- df
        proposal_list[[ct]]$state[i] <- 1
        ct <- ct + 1
        proposal_list[[ct]] <- df
        proposal_list[[ct]]$state[i] <- 2
      }
    }
    if (df$state[i] == 1) {
      ct <- ct + 1
      proposal_list[[ct]] <- df
      proposal_list[[ct]]$state[i] <- 2
    }
  }
  proposal_list
}



run_df_as_model_glmm <- function(df, glmm_dat_df, ranef) {
  
  formula <- formula_from_df_glmm(df, ranef)
  
  fit <- lme4::glmer(formula, data = glmm_dat_df, family = "binomial", control = glmerControl(optimizer = "bobyqa"))
  
  if (sum(df$state) > 1) {
    vif <- car::vif(fit)
    
    if (max(vif) > vif_threshold) return(NA)
  }
  
  return(fit)
}


initialize_modsel_df <- function(covars_base) {
  current_mod <- data.frame(
    param = covars_base,
    state = 0,
    scale = NA
  )
}

get_AIC <- function(fit) {
  suppressWarnings(if (length(fit) == 1 && is.na(fit)) return(Inf))
  
  AIC(fit)
}

do_modsel <- function(datlist, ranef = "(1|Year)", scales, ncores) {
  
  cl <- makeCluster(ncores)
  
  cap <- clusterEvalQ(cl, source("code/modsel_helper.R"))
  
  covars_base <- c("NDVI", "NDWI", "Elevation", "Slope",
              "Percent_Deciduous", "Percent_Evergreen", "Percent_Mixed",
              "Percent_Impervious")
  combos <- expand.grid(covars_base, c("_25", "_250"))
  covars <- c("Longitude", "Latitude", paste0(combos[, 1], combos[, 2]))
  covars_base  <- c("Longitude", "Latitude", covars_base)
  
  dat <- datlist[[1]]
  
  current_mod <- initialize_modsel_df(covars_base)
  
  fit_null <- run_df_as_model_glmm(df = current_mod, dat, ranef)
  # browser()
  current_AIC <- get_AIC(fit_null)
  done <- F
  ct <- 0
  
  while (!done) {
    ct <- ct + 1
    
    writeLines(paste0("Iter ", ct,  "\t--AIC: ", round(current_AIC, 2)))
    
    proposal_list <- make_proposals_glmm(current_mod, scales)

    reslist <- parLapply(cl, proposal_list, 
                         run_df_as_model_glmm, 
                         glmm_dat_df = dat,  
                         ranef = ranef)
    
    res_AIC <- reslist %>% lapply(function(x) get_AIC(x)) %>% unlist()
    
    if (all(res_AIC >= current_AIC)) {
      done <- TRUE
    } else {
      current_AIC <- min(res_AIC)
      current_mod <- proposal_list[[which(res_AIC == current_AIC)]]
    }
  }
  
  stopCluster(cl)
  
  #### Re-fit the best model ####
  
  fit_final <- run_df_as_model_glmm(current_mod, dat, ranef = ranef)
  
  return(list(mod = fit_final, 
              AIC = get_AIC(fit_final),
              best_mod_df = current_mod, 
              ranef_short = paste(ranef, collapse = ", "),
              datlist = datlist))
}
