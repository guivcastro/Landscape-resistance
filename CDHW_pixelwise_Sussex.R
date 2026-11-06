###### Pixel-wise Compound Drought-Heatwave (CDHW) Exposure - Sussex, 2022
# the cumulative daily product of standardised maximum temperature anomalies and drought 
# severity during the Summer of 2022 at 1 km resolution

library(terra)
library(sf)
library(dplyr)
library(stringr)
library(SPEI)
library(lubridate)

# working directory
setwd("C:/Users/guilh/OneDrive/Documentos/PhD/Data/Climate/EWE")

# Study area
sussex_shp  <- st_read("Sussex_boundaries.shp")
sussex_vect <- vect(sussex_shp)

# Inputs
rf_dir <- "C:/Users/guilh/OneDrive/Documentos/PhD/Data/Climate/EWE/raw_data/monthly_avg_rainfall"
tas_dir <- "C:/Users/guilh/OneDrive/Documentos/PhD/Data/Climate/EWE/raw_data/monthly_avg_temp"
tasmax_dir <- "C:/Users/guilh/OneDrive/Documentos/PhD/Data/Climate/EWE/raw_data/daily_max_temp"

###### Read monthly NetCDFs (rainfall & temperature) clipped to Sussex,

read_nc_stack_clipped <- function(pattern, path, sussex_vect) {

  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(paste("No files found for pattern:", pattern))

  # Sort files by the year found in the filename, so layer order is
  # guaranteed to be chronological before combining into one stack.
  years <- as.numeric(str_extract(basename(files), "\\d{4}"))
  files <- files[order(years)]
  years <- sort(years)
  stack_list <- list()
  index_list <- list()

  for (i in seq_along(files)) {

    r <- rast(files[i])
    crs(r) <- "EPSG:27700"
    r <- crop(r, sussex_vect)
    r <- mask(r, sussex_vect)
    n <- nlyr(r)
    names(r) <- paste0(years[i], "_", sprintf("%02d", seq_len(n)))
    stack_list[[i]] <- r
    index_list[[i]] <- data.frame(year = years[i], month = seq_len(n))
  }

  full_stack <- do.call(c, stack_list)
  index_df   <- bind_rows(index_list)
  list(stack = full_stack, index = index_df)
}

rain_res <- read_nc_stack_clipped("rainfall.*\\.nc$", rf_dir, sussex_vect)
temp_res <- read_nc_stack_clipped("tas.*\\.nc$",       tas_dir, sussex_vect)

rain_stack <- rain_res$stack
temp_stack <- temp_res$stack

# Check: rainfall and temperature must be aligned month-for-month
stopifnot(identical(rain_res$index, temp_res$index))
ym_index <- rain_res$index

# Convert kelvin to Celsius if necessary
temp_max_val <- max(global(temp_stack, "max", na.rm = TRUE)[, 1], na.rm = TRUE)
if (temp_max_val > 100) {
  temp_stack <- temp_stack - 273.15
}

###### Build a latitude raster (needed for pixel-wise Thornthwaite PET)
template <- temp_stack[[1]]

# Get every cell centre
xy <- crds(template, df = TRUE, na.rm = FALSE)
pts    <- vect(xy, geom = c("x", "y"), crs = crs(template))
pts_ll <- project(pts, "EPSG:4326")
lat_vals <- crds(pts_ll)[, 2]
lat_r <- rast(template)
values(lat_r) <- lat_vals
names(lat_r) <- "lat"

###### Pixel-wise PET (Thornthwaite)
n_months <- nlyr(temp_stack)
temp_lat_stack <- c(temp_stack, lat_r)
pet_stack <- app(temp_lat_stack, fun = function(v) {

  temp_ts <- v[1:n_months]
  lat_val <- v[n_months + 1]

  if (any(is.na(temp_ts)) || is.na(lat_val)) {
    return(rep(NA_real_, n_months))
  }

  as.numeric(thornthwaite(temp_ts, lat = lat_val))
})

names(pet_stack) <- names(temp_stack)

###### Pixel-wise water balance and SPEI-3
wb_stack <- rain_stack - pet_stack
names(wb_stack) <- names(temp_stack)
start_year  <- ym_index$year[1]
start_month <- ym_index$month[1]
spei_stack <- app(wb_stack, fun = function(v) {

  if (all(is.na(v))) return(rep(NA_real_, length(v)))

  wb_ts <- ts(v, start = c(start_year, start_month), frequency = 12)

  fit <- tryCatch(
    spei(wb_ts, scale = 3,
         ref.start = c(1991, 1),
         ref.end   = c(2020, 12))$fitted,
    error = function(e) rep(NA_real_, length(v))
  )

  as.numeric(fit)
})

names(spei_stack) <- names(temp_stack)

###### Save 2022 monthly SPEI-3 rasters
idx_2022 <- which(ym_index$year == 2022)
spei_2022 <- spei_stack[[idx_2022]]
writeRaster(
  spei_2022,
  "Sussex_SPEI3_2022_monthly.tif",
  overwrite = TRUE
)

###### Read daily Tmax NetCDFs (1991-2024)
tasmax_files <- list.files(
 tasmax_dir,
 pattern = "tasmax.*\\.nc$",
 full.names = TRUE
)

if (length(tasmax_files) == 0)
 stop("No daily Tmax files found")

tmax_layers <- list()
tmax_dates  <- list()

for (f in sort(tasmax_files)) {
 
 r <- rast(f)
 
 # Assert CRS
 crs(r) <- "EPSG:27700"
 
 r <- crop(r, sussex_vect)
 r <- mask(r, sussex_vect)
 dates <- as.Date(time(r))
 tmax_layers[[length(tmax_layers) + 1]] <- r
 tmax_dates[[length(tmax_dates) + 1]]   <- dates
}

tmax_stack <- do.call(c, tmax_layers)
tmax_dates <- do.call(c, tmax_dates)

# Order chronologically
ord <- order(tmax_dates)
tmax_stack <- tmax_stack[[ord]]
tmax_dates <- tmax_dates[ord]
names(tmax_stack) <- as.character(tmax_dates)

###### Calculate long-term SUMMER (1991-2024) per-pixel mean and sd
summer_all <- month(tmax_dates) %in% c(6, 7, 8)
tmax_summer_stack <- tmax_stack[[summer_all]]
tmax_mean_r <- app(tmax_summer_stack, fun = mean, na.rm = TRUE)
tmax_sd_r <- app(tmax_summer_stack, fun = sd, na.rm = TRUE)

###### Extract Summer 2022 (June-August)
summer2022 <- year(tmax_dates) == 2022 &
 month(tmax_dates) %in% c(6, 7, 8)
tmax_2022_stack <- tmax_stack[[summer2022]]
tmax_2022_dates <- tmax_dates[summer2022]

###### Standardise Summer 2022 Tmax using the long-term baseline
tmax_std_stack <- (tmax_2022_stack - tmax_mean_r) / tmax_sd_r
names(tmax_std_stack) <- names(tmax_2022_stack)

###### Match each summer day to its monthly SPEI-3 value
months_2022 <- month(tmax_2022_dates)
spei_daily_stack <- rast(
 lapply(months_2022, function(m) spei_2022[[m]])
)

names(spei_daily_stack) <- names(tmax_std_stack)

###### Calculate summer CDHW exposure
daily_severity_stack <- tmax_std_stack * (-1 * spei_daily_stack)
CDHW_Exposure_2022 <- app(
 daily_severity_stack,
 fun = function(v) {
  if (all(is.na(v))) return(NA_real_)
  sum(v, na.rm = TRUE)
 }
)

###### Mask and export
CDHW_Exposure_2022 <- mask(CDHW_Exposure_2022, sussex_vect)
names(CDHW_Exposure_2022) <- "CDHW_Exposure_2022"

writeRaster(
 CDHW_Exposure_2022,
 "CDHW_Exposure_summer_2022_1km.tif",
 overwrite = TRUE
)

# Plot
plot(CDHW_Exposure_2022, main = "CDHW Exposure Sussex Summer 2022")

