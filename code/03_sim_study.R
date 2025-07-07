library(terra)
library(tidyterra)
library(lme4)
set.seed(15209)

# Given the predicted occ/det probs, how much power do we have to detect
# differences in observation rate using a GLMM?

pts <- vect("data/SM_sampling_locs_GLMM.shp")

nsim <- 5000
nwindows <- 8 # 56 days of sampling, approx. 2 months
pval <- numeric(nsim)


pb <- progress::progress_bar$new(total = nsim)
this_pts <- pts %>% as.data.frame()
# Simulation if our predictions are accurate
for (i in 1:nsim) {
  pb$tick()
  
  nsuccesses <- rbinom(nrow(this_pts), nwindows, this_pts$pred_rate)
  nfailures <- nwindows - nsuccesses
  
  this_fit <- glm(cbind(nsuccesses, nfailures) ~ this_pts$type, family = "binomial")
  pval[i] <- summary(this_fit)$coefficients[2, 4]
}
power <- mean(pval < 0.05)
false_negative <- 1 - mean(pval < 0.05)


pb <- progress::progress_bar$new(total = nsim)
for (i in 1:nsim) {
  pb$tick()
  
  nsuccesses <- rbinom(nrow(pts), nwindows, sample(pts$pred_rate, replace = TRUE))
  nfailures <- nwindows - nsuccesses
  
  this_fit <- glm(cbind(nsuccesses, nfailures) ~ pts$type, family = "binomial")
  pval[i] <- summary(this_fit)$coefficients[2, 4]
}
false_positive <- mean(pval < 0.05)



