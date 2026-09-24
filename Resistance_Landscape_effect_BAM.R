# Loading packages
install.packages("mirai")
install.packages("gratia")
install.packages("terra")
install.packages("mgcv")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("gstat")
install.packages("sp")
install.packages("spdep")
library(mgcv)
library(mirai)
library(gratia)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gstat)
library(sp)
library(spdep)

# Setting working directory
setwd("C:/Users/guilh/OneDrive/Documentos/R/Complexity")

######### Loading data
# Response variable
ndvi <- rast("NDVI_resistance_2022_high_pd.tif")
names(ndvi) <- "ndvi"

# Pixel-level SE of log(Rt), propagated from ARIMA prediction
ndvi_se <- rast("NDVI_log_resistance_SE_2022_combined.tif")
names(ndvi_se) <- "ndvi_se"

# Global PC1
pc_global <- rast("PC1.tif")
names(pc_global) <- "pc"

# Moving window predictors
predictor_paths <- list(
 pc1_m9 = "PC1_mean_13x13_all.tif",
 pc1_s9 = "PC1_std_13x13_all.tif"
)
predictor_stack <- rast(lapply(predictor_paths, rast))
names(predictor_stack) <- names(predictor_paths)

# CDHW exposure
exposure <- rast("CDHW_Exposure_Summer_2022_1km_high_pd.tif")
names(exposure) <- "exposure"
NAflag(exposure) <- -9999
exposure[exposure == -9999] <- NA

# Land cover map
lcm <- rast("LCM_2023_30m.tif")
names(lcm) <- "lcm"

######### Preparing patch size data
# Calculate the area of each habitat polygon and convert it to hectares.
habitats_sf <- st_read("Land_Cover_Map_2023_parcels.shp")
habitats_sf$patch_size <- as.numeric(st_area(habitats_sf)) / 10000

template <- rast(ndvi)
patch_rast <- rasterize(vect(habitats_sf), template, field = "patch_size")
names(patch_rast) <- "patch_size"

######### Align rasters
pc_global       <- project(pc_global, ndvi, method = "bilinear")
predictor_stack <- project(predictor_stack, ndvi, method = "bilinear")
patch_rast      <- project(patch_rast, ndvi, method = "bilinear")
lcm             <- project(lcm, ndvi, method = "near")

r_stack <- c(ndvi, pc_global, predictor_stack, patch_rast, lcm)
print(names(r_stack))
stopifnot("lcm" %in% names(r_stack))

######### Dataframe conversion
df <- as.data.frame(
 r_stack,
 xy = TRUE,
 na.rm = TRUE
)

df[df == -9999] <- NA

coords <- df[, c("x", "y")]

######### Assign each 30 m pixel the exposure value and ID of its 1 km cell
df$exposure <- terra::extract(exposure, coords)[, 2]
df$exp_id   <- terra::cellFromXY(exposure, coords)

######### Land cover reclassification
reclass_lcm <- function(lc) {
 lc[lc %in% c(4, 5, 6, 7)] <- 10   #grassland
 lc[lc %in% c(9)] <- 11            #heathland
 return(lc)
}
df$lcm <- reclass_lcm(df$lcm)

######### Filtering
# Keeping only pixels with complete data and the five habitat classes of interest.
# NDVI resistance is restricted to positive values because as it will be log-transformed
dat <- df %>%
 filter(
  complete.cases(ndvi, pc, pc1_m9, pc1_s9, exposure, patch_size, lcm),
  ndvi > 0,
  lcm %in% c(1, 2, 3, 10, 11)
 )

dat$lcm    <- factor(dat$lcm, levels = c(1, 2, 3, 10, 11))

######### Transform response and patch size
#Log transformation reduces skewness
dat$log_ndvi  <- log(dat$ndvi + 1e-4)
dat$log_patch <- log1p(dat$patch_size)

######### Define GAM models
# Model 1: Ecological predictors only
# Model 2: Ecological predictors + spatial GP
# Model 3: Ecological predictors + spatial GP + CDHW exposure

# Model 1 (baseline model)
form_nospatial <- log_ndvi ~ 0 + lcm +
 s(pc, k = 10) +
 s(pc1_m9, by = lcm, k = 10) +
 s(pc1_s9, by = lcm, k = 10) +
 s(log_patch, k = 30)

# Model 2
form_spatial <- log_ndvi ~ 0 + lcm +
 s(pc, k = 10) +
 s(pc1_m9, by = lcm, k = 10) +
 s(pc1_s9, by = lcm, k = 10) +
 s(log_patch, k = 30) +
 s(x, y,
   bs = "gp",
   m = -3,
   k = 5000,
   xt = list(max.knots = 20000, seed = 42))

# Model 3
form_final <- log_ndvi ~ 0 + lcm +
 s(pc, k = 10) +
 s(pc1_m9, by = lcm, k = 10) +
 s(pc1_s9, by = lcm, k = 10) +
 s(log_patch, k = 30) +
 s(exposure, k = 30) +
 s(x, y,
   bs = "gp",
   m = -3,
   k = 5000,
   xt = list(max.knots = 20000, seed = 42))

######### Fit GAM models
# Model 1
fit_nospatial <- bam(
 form_nospatial,
 data = dat,
 method = "fREML",
 discrete = TRUE,
 nthreads = c(4, 1)
)

# Model 2
fit_spatial <- bam(
 form_spatial,
 data = dat,
 method = "fREML",
 discrete = TRUE,
 nthreads = c(4, 1)
)

# Model 3
fit_final <- bam(
 form_final,
 data = dat,
 method = "fREML",
 discrete = TRUE,
 nthreads = c(4, 1)
)


######### Extract Gaussian-process parameters from Model 2 & 3 to be compared
# Model 3
gp_id <- which(
 sapply(
  fit_final$smooth,
  inherits,
  "gp.smooth"
 )
)

gp_def <- fit_final$smooth[[gp_id]]$gp.defn

cat("Gaussian-process parameters: Final model\n")
print(gp_def)

# Model 2
gp_id_spatial <- which(
 sapply(
  fit_spatial$smooth,
  inherits,
  "gp.smooth"
 )
)

gp_def_spatial <-
 fit_spatial$smooth[[gp_id_spatial]]$gp.defn

cat("Gaussian-process parameters: Spatial-only model\n")
print(gp_def_spatial)


######### Model comparison
AIC(fit_nospatial, fit_spatial, fit_final)

######### Basis-dimension (k) sensitivity analysis
# Aim: Test whether the model results depend strongly on the chosen basis dimensions.
#
# In a GAM, k controls the maximum flexibility available to a smooth. It is not the effective degrees of freedom.
#
# We fit several models with different k values and use k.check() to assess whether the basis dimensions are sufficiently large.
#
# If increasing k does not substantially change the model, the chosen basis dimensions are likely adequate.

k_grid <- list(
 list(pc = 4, pc1 = 5, patch = 20, exposure = 10),
 list(pc = 6, pc1 = 6, patch = 25, exposure = 15),
 list(pc = 8, pc1 = 8, patch = 30, exposure = 20),
 list(pc = 10, pc1 = 10, patch = 30, exposure = 30)
)

# Storage for fitted models
fit_list <- vector(
 "list",
 length(k_grid)
)

names(fit_list) <-
 paste0("k", seq_along(k_grid))


# Storage for results
results <- data.frame(
 Model = character(),
 k_pc = integer(),
 k_pc1 = integer(),
 k_patch = integer(),
 k_exp = integer(),
 AIC = numeric(),
 stringsAsFactors = FALSE
)

# Fit sensitivity models
for (i in seq_along(k_grid)) {
 
 pars <- k_grid[[i]]
 
 form_tmp <- as.formula(
  paste0(
   "log_ndvi ~ 0 + lcm + ",
   "s(pc, k=", pars$pc, ") + ",
   "s(pc1_m9, by=lcm, k=", pars$pc1, ") + ",
   "s(pc1_s9, by=lcm, k=", pars$pc1, ") + ",
   "s(log_patch, k=", pars$patch, ") + ",
   "s(exposure, k=", pars$exposure, ") + ",
   "s(x, y, bs='gp', m=-3, k=5000, ",
   "xt=list(max.knots=20000, seed=42))"
  )
 )
 
 fit_list[[i]] <- bam(
  form_tmp,
  data = dat,
  method = "fREML",
  discrete = TRUE,
  nthreads = c(4, 1)
 )
 
 results[i, ] <- data.frame(
  Model = names(fit_list)[i],
  k_pc = pars$pc,
  k_pc1 = pars$pc1,
  k_patch = pars$patch,
  k_exp = pars$exposure
 )
 
 cat("Model", names(fit_list)[i], "\n")
 
 print(
  k.check(fit_list[[i]])
 )
}

print(results)


######### Spatial GP basis-rank sensitivity
# Aim: Evaluate whether model results remain stable as the maximum basis rank of the spatial Gaussian process increases (k = 1500, 2000, 3000, 4000, 5000)
#
# If results are similar, the spatial component is robust to the choice of basis rank.
#
# Increasing k does not force the model to become more complex, it only allows greater maximum complexity.
#
# k-index	Interpretation
# ≈ 1	No strong evidence that k is too small
# < 1	Potential evidence that k may be too small
# more negative ( <1) indicates that the basis may be inadequate

gp_grid <- c(
 1500,
 2000,
 3000,
 4000,
 5000
)

# Storage for fitted models
gp_sensitivity <- vector(
 "list",
 length(gp_grid)
)

names(gp_sensitivity) <- paste0(
 "GP_",
 gp_grid
)

# Storage for diagnostic results
gp_results <- data.frame(
 GP_k = numeric(),
 GP_EDF = numeric(),
 GP_k_index = numeric(),
 GP_p_value = numeric(),
 stringsAsFactors = FALSE
)

# Fit GP sensitivity models
for (i in seq_along(gp_grid)) {
 
 k_gp <- gp_grid[i]
 
 cat("GP k =", k_gp, "\n")
 
 # Model formula
 form_gp <- as.formula(
  paste0(
   "log_ndvi ~ 0 + lcm + ",
   "s(pc, k=10) + ",
   "s(pc1_m9, by=lcm, k=10) + ",
   "s(pc1_s9, by=lcm, k=10) + ",
   "s(log_patch, k=30) + ",
   "s(exposure, k=30) + ",
   "s(x, y, bs='gp', m=-3, k=",
   k_gp,
   ", xt=list(max.knots=20000, seed=42))"
  )
 )
 
 # Fit model
 gp_sensitivity[[i]] <- bam(
  form_gp,
  data = dat,
  method = "fREML",
  discrete = TRUE,
  nthreads = c(4, 1)
 )
 
 # Extract GP EDF
 gp_table <- summary(
  gp_sensitivity[[i]]
 )$s.table
 gp_row <- grep(
  "^s\\(x,y\\)",
  rownames(gp_table)
 )
 gp_edf <- gp_table[
  gp_row,
  "edf"
 ]
 
 # k.check() diagnostics
 k_tab <- k.check(
  gp_sensitivity[[i]]
 )
 gp_k_row <- grep(
  "^s\\(x,y\\)",
  rownames(k_tab)
 )
 
 if (length(gp_k_row) == 1) {
  gp_k_index <- k_tab[
   gp_k_row,
   "k-index"
  ]
  gp_p_value <- k_tab[
   gp_k_row,
   "p-value"
  ]
 } else {
  gp_k_index <- NA
  gp_p_value <- NA
 }
 
 # Store results
 gp_results <- rbind(
  gp_results,
  data.frame(
   GP_k = k_gp,
   GP_EDF = gp_edf,
   GP_k_index = gp_k_index,
   GP_p_value = gp_p_value
  )
 )
 
 # Print basic output
 cat(
  "GP EDF:",
  round(gp_edf, 2),
  "\n"
 )
 
 cat(
  "GP k-index:",
  round(gp_k_index, 3),
  "\n"
 )
 
 cat(
  "GP p-value:",
  round(gp_p_value, 4),
  "\n"
 )
 
}

# Final supplementary-results table
cat("GP basis-rank sensitivity summary\n")

print(
 gp_results,
 row.names = FALSE
)

######### Calculate model residuals
# Residuals represent the part of NDVI resistance that the model does not explain.

dat$resid_nospatial <- residuals(fit_nospatial, type = "deviance")
dat$resid_spatial   <- residuals(fit_spatial, type = "deviance")
dat$resid_final     <- residuals(fit_final, type = "deviance")

######### Moran's I - residual spatial autocorrelation
# Aim: Test amongst GAMs whether spatial autocorrelation remains in the residuals.
# The same 50,000 pixels are used for all models so that differences in Moran's I are caused by the models rather than different samples.

set.seed(42)

samp_id <- sample(seq_len(nrow(dat)), 50000)

samp <- dat[samp_id, ]

# Convert to spatial points
coordinates(samp) <- ~x + y

coords_samp <- coordinates(samp)

# 8-nearest-neighbour spatial weights
nb <- knearneigh(coords_samp, k = 8)

lw <- nb2listw(
 knn2nb(nb),
 style = "W"
)

# Moran's I: no spatial term
cat("Moran's I: No spatial term\n")
print(
 moran.test(
  samp$resid_nospatial,
  lw
 )
)

# Moran's I: spatial GP
cat("Moran's I: Spatial GP\n")
print(
 moran.test(
  samp$resid_spatial,
  lw
 )
)

# Moran's I: final model
cat("Moran's I: Final model\n")
print(
 moran.test(
  samp$resid_final,
  lw
 )
)

# Residual variograms
res_sample <- data.frame(
 x = samp$x,
 y = samp$y,
 resid_nospatial = samp$resid_nospatial,
 resid_spatial = samp$resid_spatial,
 resid_final = samp$resid_final
)

# Model 1
v0 <- variogram(
 resid_nospatial ~ 1,
 locations = ~x + y,
 data = res_sample,
 cutoff = 5000,
 width = 100
)

# Model 2
v1 <- variogram(
 resid_spatial ~ 1,
 locations = ~x + y,
 data = res_sample,
 cutoff = 5000,
 width = 100
)

# Model 3
v2 <- variogram(
 resid_final ~ 1,
 locations = ~x + y,
 data = res_sample,
 cutoff = 5000,
 width = 100
)

# Combine for plotting
v0$model <- "No spatial term"
v1$model <- "Spatial GP"
v2$model <- "Final model"

variogram_df <- rbind(
 v0[, c("dist", "gamma", "model")],
 v1[, c("dist", "gamma", "model")],
 v2[, c("dist", "gamma", "model")]
)

# Plot residual variograms
ggplot(
 variogram_df,
 aes(x = dist, y = gamma, linetype = model)
) +
 geom_line(linewidth = 1) +
 labs(
  x = "Distance (m)",
  y = "Residual semivariance",
  linetype = NULL
 ) +
 theme_bw() +
 theme(
  panel.grid = element_blank()
 )

######### Final model diagnostics
summary(fit_final)
gam.check(fit_final)
k.check(fit_final)

######### Concurvity
# Aim: Test whether predictors contain overlapping information (50,000 pixels only)
set.seed(42)

dat_con <- dat[sample(nrow(dat), 50000), ]

form_con <- log_ndvi ~ 0 + lcm +
 s(pc, k = 10) +
 s(pc1_m9, by = lcm, k = 10) +
 s(pc1_s9, by = lcm, k = 10) +
 s(log_patch, k = 30) +
 s(exposure, k = 30) +
 s(x, y, bs = "gp", m = -3, k = 1000,
   xt = list(max.knots = 20000, seed = 42))

fit_con <- bam(
 form_con,
 data = dat_con,
 method = "fREML",
 discrete = TRUE,
 nthreads = c(4, 1)
)

conc <- concurvity(fit_con, full = FALSE)

conc

######### Plotting habitat-specific smooth terms altogether
lcm_labels <- c(
 "1"  = "Broadleaved woodland",
 "2"  = "Coniferous woodland",
 "3"  = "Arable",
 "10" = "Grassland",
 "11" = "Heathland"
)
habitats <- names(lcm_labels)

extract_smooth <- function(var, habitat) {
 term <- paste0("s(", var, "):lcm", habitat)
 
 newdat <- data.frame(
  pc = mean(dat$pc, na.rm = TRUE),
  
  pc1_m9 = if (var == "pc1_m9")
   seq(min(dat$pc1_m9, na.rm = TRUE), max(dat$pc1_m9, na.rm = TRUE), length.out = 200)
  else mean(dat$pc1_m9, na.rm = TRUE),
  
  pc1_s9 = if (var == "pc1_s9")
   seq(min(dat$pc1_s9, na.rm = TRUE), max(dat$pc1_s9, na.rm = TRUE), length.out = 200)
  else mean(dat$pc1_s9, na.rm = TRUE),
  
  log_patch = mean(dat$log_patch, na.rm = TRUE),
  x = mean(dat$x, na.rm = TRUE),
  y = mean(dat$y, na.rm = TRUE),
  
  exposure = mean(dat$exposure, na.rm = TRUE),
  lcm = factor(habitat, levels = levels(dat$lcm))
 )
 
 pred <- predict(
  fit_final,
  newdata = newdat,
  type = "terms",
  se.fit = TRUE
 )
 
 data.frame(
  habitat = habitat,
  x = if (var == "pc1_m9") newdat$pc1_m9 else newdat$pc1_s9,
  fit = pred$fit[, term],
  se = pred$se.fit[, term],
  type = var
 )
}

plot_df <- bind_rows(
 lapply(habitats, function(h) extract_smooth("pc1_m9", h)),
 lapply(habitats, function(h) extract_smooth("pc1_s9", h))
)

plot_df <- plot_df %>%
 mutate(
  lower = fit - 2 * se,
  upper = fit + 2 * se,
  habitat = recode(habitat, !!!lcm_labels),
  type = recode(type, "pc1_m9" = "Mean PC1", "pc1_s9" = "SD PC1")
 )

habitat_colours <- c(
 "Broadleaved woodland" = "#00734C",
 "Coniferous woodland"  = "#6FE64E",
 "Arable"               = "#8B4513",
 "Grassland"            = "#F2B300",
 "Heathland"            = "#800080"
)

x_breaks_mean <- scale_x_continuous(
 breaks = seq(floor(min(plot_df$x, na.rm = TRUE)),
              ceiling(max(plot_df$x, na.rm = TRUE)), by = 1)
)
x_breaks_std <- scale_x_continuous(
 breaks = seq(floor(min(plot_df$x, na.rm = TRUE) * 2) / 2,
              ceiling(max(plot_df$x, na.rm = TRUE) * 2) / 2, by = 0.5)
)

base_theme <- theme_bw(base_size = 13) +
 theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 20),
  strip.background = element_blank(),
  legend.position = "bottom"
 )

p_mean <- ggplot(subset(plot_df, type == "Mean PC1"),
                 aes(x = x, y = fit, colour = habitat, fill = habitat)) +
 geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
 geom_line(linewidth = 1) +
 geom_hline(yintercept = 0, linetype = 2) +
 scale_colour_manual(values = habitat_colours) +
 scale_fill_manual(values = habitat_colours) +
 guides(colour = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
 base_theme +
 labs(title = "Landscape composition (mean of PC1)", x = "PC1", y = "Log(Rt)") +
 x_breaks_mean

p_std <- ggplot(subset(plot_df, type == "SD PC1"),
                aes(x = x, y = fit, colour = habitat, fill = habitat)) +
 geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
 geom_line(linewidth = 1) +
 geom_hline(yintercept = 0, linetype = 2) +
 scale_colour_manual(values = habitat_colours) +
 scale_fill_manual(values = habitat_colours) +
 guides(colour = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
 base_theme +
 labs(title = "Landscape heterogeneity", x = "PC1", y = "Log(Rt)") +
 x_breaks_std

final_plot <- patchwork::wrap_plots(p_mean, p_std, ncol = 2)
final_plot

ggsave(
 filename = "PC1_plots.png",
 plot = final_plot,
 width = 12, height = 6, dpi = 330, units = "in", bg = "white"
)


###### Plotting habitat-specific smooth terms
library(ggplot2)
library(dplyr)
library(purrr)
library(cowplot)

# Habitat labels
lcm_labels <- c(
 "1"  = "Broadleaved woodland",
 "2"  = "Coniferous woodland",
 "3"  = "Arable",
 "10" = "Grassland",
 "11" = "Heathland"
)

habitats <- names(lcm_labels)

# Extract smooths
extract_smooth <- function(var, habitat) {
 
 term <- paste0("s(", var, "):lcm", habitat)
 
 newdat <- data.frame(
  
  pc = mean(
   dat$pc,
   na.rm = TRUE
  ),
  
  pc1_m9 = if (var == "pc1_m9") {
   
   seq(
    min(dat$pc1_m9, na.rm = TRUE),
    max(dat$pc1_m9, na.rm = TRUE),
    length.out = 200
   )
   
  } else {
   mean(
    dat$pc1_m9,
    na.rm = TRUE
   )
  },
  
  pc1_s9 = if (var == "pc1_s9") {
   
   seq(
    min(dat$pc1_s9, na.rm = TRUE),
    max(dat$pc1_s9, na.rm = TRUE),
    length.out = 200
   )
  } else {
   mean(
    dat$pc1_s9,
    na.rm = TRUE
   )
  },
  
  log_patch = mean(
   dat$log_patch,
   na.rm = TRUE
  ),
  
  x = mean(
   dat$x,
   na.rm = TRUE
  ),
  
  y = mean(
   dat$y,
   na.rm = TRUE
  ),
  
  exposure = mean(
   dat$exposure,
   na.rm = TRUE
  ),
  
  lcm = factor(
   habitat,
   levels = levels(dat$lcm)
  )
 )
 
 # Predictions
 pred <- predict(
  fit_final,
  newdata = newdat,
  type = "terms",
  se.fit = TRUE
 )
 
 data.frame(
  habitat = habitat,
  x = if (var == "pc1_m9") {
   newdat$pc1_m9
  } else {
   newdat$pc1_s9
  },
  fit = pred$fit[, term],
  se = pred$se.fit[, term],
  type = var
 )
}

# Extract all curves
plot_df <- bind_rows(
 lapply(
  habitats,
  function(h) {
   extract_smooth(
    "pc1_m9",
    h
   )
  }
 ),
 
 lapply(
  habitats,
  function(h) {
   extract_smooth(
    "pc1_s9",
    h
   )
  }
 )
)

# Confidence intervals
plot_df <- plot_df %>%
 
 mutate(
  lower = fit - 2 * se,
  upper = fit + 2 * se,
  habitat = recode(
   habitat,
   !!!lcm_labels
  ),
  
  type = recode(
   type,
   "pc1_m9" = "Mean PC1",
   "pc1_s9" = "SD PC1"
  )
 )

# Habitat colours
habitat_colours <- c(
 "Broadleaved woodland" = "#00734C",
 "Coniferous woodland" = "#6FE64E",
 "Arable" = "#8B4513",
 "Grassland" = "#F2B300",
 "Heathland" = "#800080"
)

# Base theme
base_theme <- theme_bw(
 base_size = 13
) +
 
 theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(
   colour = "black"
  ),
  # AXIS NUMBERS
  axis.text.x = element_text(
   size = 25
  ),
  axis.text.y = element_text(
   size = 25
  ),
  # AXIS TITLES
  axis.title.x = element_text(
   size = 20
  ),
  axis.title.y = element_text(
   size = 20
  ),
  # PLOT TITLE
  plot.title = element_text(
   size = 18,
   face = "bold",
   hjust = 0.5
  ),
  legend.position = "none",
  plot.margin = margin(
   8,
   8,
   8,
   8
  )
 )

# Function to create individual habitat plots
make_habitat_plot <- function(
  data,
  habitat_name,
  x_label,
  x_breaks
) {
 
 ggplot(
  data,
  aes(
   x = x,
   y = fit
  )
 ) +
 
 geom_ribbon(
  aes(
   ymin = lower,
   ymax = upper
  ),
  fill = habitat_colours[habitat_name],
  alpha = 0.18,
  colour = NA
 ) +
 
 geom_line(
  colour = habitat_colours[habitat_name],
  linewidth = 1.2
 ) +
 
 geom_hline(
  yintercept = 0,
  linetype = 2
 ) +
 
 scale_x_continuous(
  breaks = x_breaks,
  labels = x_breaks,
  expand = expansion(
   mult = c(
    0.02,
    0.02
   )
  )
 ) +
 
 labs(
  title = habitat_name,
  x = x_label,
  y = "Log(Rt)"
 ) +
  
 base_theme
}

# Mean PC1
mean_range <- range(
 plot_df$x[
  plot_df$type == "Mean PC1"
 ],
 na.rm = TRUE
)

mean_breaks <- pretty(
 mean_range,
 n = 5
)

# Keep only breaks inside the actual range
mean_breaks <- mean_breaks[
 mean_breaks >= mean_range[1] &
  mean_breaks <= mean_range[2]
]

# SD PC1
std_range <- range(
 plot_df$x[
  plot_df$type == "SD PC1"
 ],
 na.rm = TRUE
)

std_breaks <- pretty(
 std_range,
 n = 5
)

# Keep only breaks inside actual range
std_breaks <- std_breaks[
 std_breaks >= std_range[1] &
  std_breaks <= std_range[2]
]

# Create mean PC1 Panels
mean_panels <- lapply(
 
 habitats,
 
 function(h) {
  
  h_name <- lcm_labels[h]
  
  make_habitat_plot(
   
   data = plot_df %>%
    filter(
     habitat == h_name,
     type == "Mean PC1"
    ),
   
   habitat_name = h_name,
   
   x_label = "Mean PC1",
   
   x_breaks = mean_breaks
  )
 }
)


# Create SD PC1 panels
std_panels <- lapply(
 
 habitats,
 
 function(h) {
  
  h_name <- lcm_labels[h]
  
  make_habitat_plot(
   
   data = plot_df %>%
    filter(
     habitat == h_name,
     type == "SD PC1"
    ),
   
   habitat_name = h_name,
   
   x_label = "SD PC1",
   
   x_breaks = std_breaks
  )
 }
)

# Display each plot
mean_panels[[1]]
mean_panels[[2]]
mean_panels[[3]]
mean_panels[[4]]
mean_panels[[5]]

std_panels[[1]]
std_panels[[2]]
std_panels[[3]]
std_panels[[4]]
std_panels[[5]]

# Create output directories
dir.create(
 "habitat_plots_mean",
 showWarnings = FALSE
)

dir.create(
 "habitat_plots_sd",
 showWarnings = FALSE
)

# File names
habitat_names <- c(
 "Broadleaved_woodland",
 "Coniferous_woodland",
 "Arable",
 "Grassland",
 "Heathland"
)

# Save mean PC1 Plots
for (i in seq_along(mean_panels)) {
 
 ggsave(
  filename = file.path(
   "habitat_plots_mean",
   paste0(
    habitat_names[i],
    "_Mean_PC1.png"
   )
  ),
  
  plot = mean_panels[[i]],
  
  width = 6.5,
  height = 3,
  
  units = "in",
  
  dpi = 300
 )
}

# Save SD PC1 Plots
for (i in seq_along(std_panels)) {
 
 ggsave(
  filename = file.path(
   "habitat_plots_sd",
   paste0(
    habitat_names[i],
    "_SD_PC1.png"
   )
  ),
  
  plot = std_panels[[i]],
  
  width = 6.5,
  height = 3,
  
  units = "in",
  
  dpi = 300
 )
}


