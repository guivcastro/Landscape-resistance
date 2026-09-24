# Load libraries
library(terra)      # raster handling
library(sf)         # shapefiles
library(dplyr)      # data manipulation
library(stringr)    # string handling
library(SPEI)       # SPEI computation
library(ggplot2)
library(lubridate)

# Set working directory
setwd("C:/PATH/TO/WD")

# Load Sussex shapefile
sussex_shp <- st_read("Sussex_boundaries.shp")
sussex_vect <- vect(sussex_shp)

# Function to read HadUK-Grid NetCDF folder and compute AOI mean
read_nc_folder_clipped <- function(pattern, path = ".") {
 
 files <- list.files(path, pattern = pattern, full.names = TRUE)
 if(length(files) == 0) stop(paste("No files found for pattern:", pattern))
 
 data_list <- lapply(files, function(f) {
  
  r <- rast(f)
  
  # Clip to Sussex AOI
  r_sussex <- crop(r, sussex_vect)
  r_sussex <- mask(r_sussex, sussex_vect)
  
  monthly_mean <- global(r_sussex, "mean", na.rm = TRUE)[,1]
  
  year <- as.numeric(str_extract(basename(f), "\\d{4}"))
  
  data.frame(
   year  = rep(year, 12),
   month = 1:12,
   value = monthly_mean
  )
 })
 
 bind_rows(data_list) %>% arrange(year, month)
}

# Paths
rf_dir  <- "C:/PATH/TO/FILES/"
tas_dir <- "C:/PATH/TO/FILES/"

# Load data (AOI-averaged)
rain_sussex <- read_nc_folder_clipped("rainfall.*\\.nc$", rf_dir)
temp_sussex <- read_nc_folder_clipped("tas.*\\.nc$", tas_dir)

# Convert Kelvin → Celsius if necessary
if(max(temp_sussex$value, na.rm = TRUE) > 100) {
 temp_sussex$value <- temp_sussex$value - 273.15
}

############### Drought Event Detection #####################
# PET (Thornthwaite)
lat_centroid <- st_coordinates(st_centroid(st_geometry(sussex_shp)))[,2]

pet <- thornthwaite(
 temp_sussex$value,
 lat = mean(lat_centroid)
)

# Water balance
wb <- rain_sussex$value - pet

wb_ts <- ts(
 data = wb,
 start = c(1991, 1),
 frequency = 12
)

# SPEI-3 (1991–2020 climatology)
spei_3 <- spei(
 wb_ts,
 scale = 3,
 ref.start = c(1991, 1),
 ref.end   = c(2020, 12)
)

spei_df <- data.frame(
 year  = rain_sussex$year,
 month = rain_sussex$month,
 SPEI3 = as.numeric(spei_3$fitted)
)

threshold <- -1      # Moderate drought
min_duration <- 2    # Recommended for event studies

spei_df <- spei_df %>%
 mutate(drought_month = SPEI3 <= threshold)

# Run-length encoding
r <- rle(spei_df$drought_month)

spei_df$event_id <- rep(seq_along(r$lengths), r$lengths)
spei_df$event_id[!spei_df$drought_month] <- NA

# Compute duration
spei_df <- spei_df %>%
 group_by(event_id) %>%
 mutate(duration = ifelse(!is.na(event_id), n(), NA)) %>%
 ungroup()

# Apply minimum duration
spei_df <- spei_df %>%
 mutate(
  drought_event = !is.na(event_id) & duration >= min_duration,
  event_id = ifelse(drought_event, event_id, NA)
 )

############## Event Metrics
event_summary <- spei_df %>%
 filter(drought_event) %>%
 group_by(event_id) %>%
 summarise(
  start_year  = first(year),
  start_month = first(month),
  end_year    = last(year),
  end_month   = last(month),
  duration    = n(),
  severity    = sum(abs(SPEI3)),
  intensity   = min(SPEI3),
  .groups = "drop"
 )

# Save outputs
write.csv(spei_df, "Sussex_monthly_SPEI3_1991_2024.csv", row.names = FALSE)
write.csv(event_summary, "Sussex_drought_events.csv", row.names = FALSE)

# Histogram
hist(
 spei_df$SPEI3,
 breaks = 60,
 col = "skyblue",
 border = "white",
 xlab = "SPEI-3",
 ylab = "Frequency"
)


# Time-series plot
spei_df$Date <- as.Date(paste(spei_df$year, spei_df$month, "01", sep = "-"))

x_range <- range(spei_df$Date)
y_range <- range(spei_df$SPEI3)

ggplot(spei_df, aes(x = Date, y = SPEI3, fill = SPEI3)) +
 geom_col(width = 25) +
 scale_fill_gradient2(
  low = "#7A4A2E",
  mid = "gray90",
  high = "blue",
  midpoint = 0,
  limits = c(-3, 3)
 ) +
 geom_hline(
  yintercept = c(-1, 1),
  color = "black",
  linetype = "dashed",
  linewidth = 0.2
 ) +
 scale_x_date(limits = x_range) +
 scale_y_continuous(limits = y_range) +
 labs(
  x = "Year",
  y = "SPEI-3",
  fill = "SPEI-3"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.3)
 )

ggsave(
 filename = "SPEI3_Sussex.tif",
 plot = last_plot(),
 device = "tiff",
 dpi = 330,
 width = 18,
 height = 10,
 units = "cm"
)


######## Heatwaves ################
# Inspect one .nc file
f <- list.files(
 "C:/PATH/TO/FILES/",
 pattern = "tasmax.*\\.nc$",
 full.names = TRUE
)[1]

r <- rast(f)
r
nlyr(r)
time(r)

# Function to read .nc files and all layers (days) within them
read_daily_tmax <- function(path, pattern = "tasmax.*\\.nc$") {
 
 files <- list.files(path, pattern = pattern, full.names = TRUE)
 if (length(files) == 0) stop("No NetCDF files found")
 files <- sort(files)
 
 out <- lapply(files, function(f) {
  
  r <- rast(f)
  
  # project AOI once per file
  sussex_proj <- project(sussex_vect, crs(r))
  
  # crop + mask
  r <- crop(r, sussex_proj)
  r <- mask(r, sussex_proj)
  
  # spatial mean for EACH layer (each day)
  daily_mean <- global(r, "mean", na.rm = TRUE)[, 1]
  
  # extract dates from NetCDF time dimension
  dates <- as.Date(time(r))
  
  # safety check
  if (length(dates) != length(daily_mean)) {
   stop("Date / layer mismatch in: ", basename(f))
  }
  
  data.frame(
   date = dates,
   Tmax = daily_mean
  )
 })
 
 bind_rows(out) %>% arrange(date)
}

#Load all daily data
tmax_df <- read_daily_tmax(
 "C:/PATH/TO/FILES/"
)

# Check
nrow(tmax_df)
range(tmax_df$date)
head(tmax_df, 10)
tail(tmax_df, 10)

#Plot daily maximum temperature
plot(
 tmax_df$date,
 tmax_df$Tmax,
 type = "l",
 col = "firebrick",
 xlab = "Date",
 ylab = "Daily mean Tmax (°C)",
)

# Add time variables
tmax_df <- tmax_df %>%
 mutate(
  year = year(date),
  doy  = yday(date)
 )

# Build daily climatological baseline (1991–2020)
ref_period <- tmax_df %>%
 filter(year >= 1991, year <= 2020)

# Compute rolling-window percentiles
climatology <- lapply(1:366, function(d) {
 
 window <- ((d - 15):(d + 15) - 1) %% 366 + 1
 
 data.frame(
  doy = d,
  t90 = quantile(
   ref_period$Tmax[ref_period$doy %in% window],
   probs = 0.9,
   na.rm = TRUE
  )
 )
}) %>%
 bind_rows()

# Attach threshold to every day
tmax_df <- tmax_df %>%
 left_join(climatology, by = "doy") %>%
 mutate(heatwave_day = Tmax > t90)

table(tmax_df$heatwave_day)

# Bridge gaps of <= 2 days so events are independent only if separated by >= 3 days
hw <- tmax_df$heatwave_day
r <- rle(hw)

lengths <- r$lengths
values  <- r$values

for(i in seq_along(values)) {
 
 if(!values[i] && lengths[i] <= 2) {
  
  if(i > 1 && i < length(values)) {
   if(values[i-1] && values[i+1]) {
    values[i] <- TRUE
   }
  }
 }
}

hw_bridged <- inverse.rle(list(lengths = lengths, values = values))

tmax_df$heatwave_day <- hw_bridged

# Heatwave detection
# Detect heatwave events (≥ 3 consecutive days)
r <- rle(tmax_df$heatwave_day)

ends   <- cumsum(r$lengths)
starts <- ends - r$lengths + 1

events <- which(r$values & r$lengths >= 3)

# Build heatwave event table
heatwave_events <- lapply(seq_along(events), function(i) {
 
 idx <- events[i]
 s <- starts[idx]
 e <- ends[idx]
 
 df <- tmax_df[s:e, ]
 
 # keep only positive exceedances
 exceed <- pmax(df$Tmax - df$t90, 0)
 
 data.frame(
  event_id   = i,
  start_date = min(df$date),
  end_date   = max(df$date),
  duration   = nrow(df),
  Tmax_mean  = mean(df$Tmax),
  Tmax_max   = max(df$Tmax),
  intensity  = mean(exceed),
  severity   = sum(exceed)
 )
}) %>% bind_rows()

# Summary of heatwave events
heatwave_events
summary(heatwave_events$duration)

# Plotting heatwave timeseries
plot(
 tmax_df$date,
 tmax_df$Tmax,
 type = "l",
 col = "grey70",
 xlab = "Years",
 ylab = "Daily maximum temperature (°C)",
)

hw_days <- tmax_df %>% filter(heatwave_day)

points(
 hw_days$date,
 hw_days$Tmax,
 col = "red",
 pch = 16
)

# Plotting heatwave durations vs intensity
# Duration vs Intensity with custom axes covering full data extent
xmax <- max(heatwave_events$duration, na.rm = TRUE)
ymax <- max(heatwave_events$intensity, na.rm = TRUE)

plot(
 heatwave_events$duration,
 heatwave_events$intensity,
 pch = 16,
 col = "darkorange",
 xlab = "Duration (days)",
 ylab = "Mean intensity (°C above threshold)",
 cex.lab = 1.5,   # axis label size
 cex.axis = 1.3,  # tick label size
 bty = "n",       # remove box
 xaxt = "n",      # suppress default x-axis
 yaxt = "n",      # suppress default y-axis
 xlim = c(0, xmax),
 ylim = c(0, ymax)
)

# Custom X-axis: every day
axis(1, at = 0:xmax, cex.axis = 1.3)

# Custom Y-axis: integer sequence from 0 to max
axis(2, at = seq(0, ceiling(ymax), by = 1), cex.axis = 1.3)

# Histogram - Number of heatwaves
# Prepare annual counts
annual_hw <- heatwave_events %>%
 mutate(year = year(start_date)) %>%
 count(year)

# Plot histogram
plot(
 annual_hw$year,
 annual_hw$n,
 type = "h",                   # vertical bars
 lwd = 10,                     # make bars thicker / wider
 col = "#964B00",            # bar color
 xlab = "Year",
 ylab = "Number of heatwaves",
 cex.lab = 1.5,                # axis label size
 cex.axis = 1.3,               # tick label size
 bty = "n",                     # remove box
 xlim = c(min(annual_hw$year), max(annual_hw$year)),
 ylim = c(0, max(annual_hw$n))
)


# Plotting heatwave durations vs severity
# Duration vs Severity with full axis coverage
xmax <- max(heatwave_events$duration, na.rm = TRUE)
ymax <- max(heatwave_events$severity, na.rm = TRUE)

plot(
 heatwave_events$duration,
 heatwave_events$severity,
 pch = 16,
 col = "darkred",         
 xlab = "Duration (days)",
 ylab = "Cumulative severity (°C-days)",
 cex.lab = 1.5,            
 cex.axis = 1.3,           
 bty = "n",                
 xaxt = "n",               
 yaxt = "n",               
 xlim = c(0, xmax),
 ylim = c(0, ymax)
)

# Custom X-axis: every day
axis(1, at = 0:xmax, cex.axis = 1.3)

# Custom Y-axis: integers from 0 to max severity
axis(2, at = seq(0, ceiling(ymax), by = 5), cex.axis = 1.3)


############ Compound drought–heatwave (CDHW) #############

# Align SPEI-3 to daily heatwave data
## Add month + year to daily Tmax
tmax_df <- tmax_df %>%
 mutate(
  year  = year(date),
  month = month(date)
 )

## Join SPEI-3
tmax_df <- tmax_df %>%
 left_join(
  spei_df %>% select(year, month, SPEI3, drought_event),
  by = c("year", "month")
 )

# Check
table(tmax_df$drought_event)
summary(tmax_df$SPEI3)

# Identify comdrought_event.x# Identify compound days
tmax_df <- tmax_df %>%
 mutate(
  CDHW_day = heatwave_day & drought_event
 )

table(tmax_df$CDHW_day)

# Identify Compound Drought–Heatwave events
tmax_df$hw_event_id <- NA_integer_

for (i in seq_len(nrow(heatwave_events))) {
 idx <- tmax_df$date >= heatwave_events$start_date[i] &
  tmax_df$date <= heatwave_events$end_date[i]
 tmax_df$hw_event_id[idx] <- heatwave_events$event_id[i]
}

# Summarise CDHW events
CDHW_events <- tmax_df %>%
 filter(!is.na(hw_event_id)) %>%
 group_by(hw_event_id) %>%
 summarise(
  start_date = min(date),
  end_date   = max(date),
  duration   = n(),
  drought_days = sum(drought_event),
  compound_days = sum(CDHW_day),
  Tmax_mean  = mean(Tmax),
  Tmax_max   = max(Tmax),
  severity_heat = sum(Tmax - t90),
  mean_SPEI3 = mean(SPEI3),
  .groups = "drop"
 ) %>%
 filter(compound_days > 0)

# Check
CDHW_events
summary(CDHW_events$compound_days)

# Define event-based columns
tmax_df <- tmax_df %>%
 mutate(
  heatwave_event_day = !is.na(hw_event_id),        # only ≥3-day heatwave days
  CDHW_event_day = heatwave_event_day & drought_event   # only compound days within heatwaves
 )

###### Frequency of CDHW
# Add year column based on event start
CDHW_events <- CDHW_events %>%
 mutate(year = year(start_date))

# Count events per year
annual_CDHw_freq <- CDHW_events %>%
 count(year, name = "n_CDHw_events")

annual_CDHw_freq

mean_annual_CDHw <- annual_CDHw_freq %>%
 summarise(mean_events_per_year = mean(n_CDHw_events)) %>%
 pull(mean_events_per_year)

mean_annual_CDHw

write.csv(
 CDHW_events,
 "CDHW_events.csv",
 row.names = FALSE
)

plot(
 annual_CDHw_freq$year,
 annual_CDHw_freq$n_CDHw_events,
 type = "h",
 lwd = 10,
 col = "#8B0000",
 xlab = "Year",
 ylab = "Number of CDHW events",
 cex.lab = 1.5,
 cex.axis = 1.1,    
 bty = "n",
 xaxt = "n",        
 xlim = c(1991, 2024),
 ylim = c(0, max(annual_CDHw_freq$n_CDHw_events))
)

axis(
 1,
 at = 1991:2024,    
 las = 1,          
 cex.axis = 0.8
)


# Plot
# Increase right margin
par(mar = c(5, 4, 4, 10))  # bottom, left, top, right

plot(
 tmax_df$date,
 tmax_df$Tmax,
 type = "l",
 col = "grey60",
 xlab = "Years",
 ylab = "Daily maximum temperature (°C)",
 bty  = "l"   
)

# heatwave events
points(
 tmax_df$date[tmax_df$heatwave_event_day],
 tmax_df$Tmax[tmax_df$heatwave_event_day],
 col = "red",
 pch = 16,
 cex = 0.7
)

# CDHW events
points(
 tmax_df$date[tmax_df$CDHW_event_day],
 tmax_df$Tmax[tmax_df$CDHW_event_day],
 col = "darkred",
 pch = 16,
 cex = 0.7
)

# Add legend
legend(
 x = max(tmax_df$date) + 50,   
 y = max(tmax_df$Tmax, na.rm = TRUE),
 legend = c("Daily max temperature", "Heatwave", "CDHW"),
 col = c("grey60", "red", "darkred"),
 pch = c(NA, 16, 16),
 lty = c(1, NA, NA),
 pt.cex = 0.7,
 bty = "n",
 cex = 0.8,
 xpd = TRUE                   
)

### Year of 2022
tmax_2022 <- tmax_df %>% filter(year(date) == 2022)

# Plot
plot(
 tmax_2022$date, tmax_2022$Tmax,
 type = "l", col = "grey60",
 xlab = "", ylab = "Daily maximum temperature (°C)",
 bty  = "l", xaxt = "n"
)

# Custom x-axis
axis.Date(side = 1, at = seq(as.Date("2022-01-01"), as.Date("2022-12-01"), by = "month"), format = "%b")
mtext("2022", side = 1, line = 2.5)

# Heatwave days
points(
 tmax_2022$date[tmax_2022$heatwave_event_day],
 tmax_2022$Tmax[tmax_2022$heatwave_event_day],
 col = "red", pch = 16, cex = 0.7
)

# CDHW days
points(
 tmax_2022$date[tmax_2022$CDHW_event_day],
 tmax_2022$Tmax[tmax_2022$CDHW_event_day],
 col = "darkred", pch = 16, cex = 0.7
)

# Add legend
legend(
 "topright",                      
 legend = c("Daily max temperature", "Heatwave", "CDHW"),
 col = c("grey60", "red", "darkred"),
 pch = c(NA, 16, 16),             
 lty = c(1, NA, NA),             
 pt.cex = 0.7,
 bty = "n",
 cex = 0.8
)


#### Exporting plots
tiff(
 filename = "CDHW_full_period.tiff",
 width = 3000,      # width in pixels
 height = 1400,     # height in pixels
 res = 330          # dpi
)

# Increase right margin
par(mar = c(5, 4, 4, 10))

plot(
 tmax_df$date,
 tmax_df$Tmax,
 type = "l",
 col = "grey60",
 xlab = "Years",
 ylab = "Daily maximum temperature (°C)",
 bty  = "l"
)

# Heatwave days
points(
 tmax_df$date[tmax_df$heatwave_event_day],
 tmax_df$Tmax[tmax_df$heatwave_event_day],
 col = "red",
 pch = 16,
 cex = 0.5
)

# CDHW days
points(
 tmax_df$date[tmax_df$CDHW_event_day],
 tmax_df$Tmax[tmax_df$CDHW_event_day],
 col = "darkred",
 pch = 16,
 cex = 0.5
)

# Legend outside plot
legend(
 x = max(tmax_df$date) + 50,
 y = max(tmax_df$Tmax, na.rm = TRUE),
 legend = c("Daily max temperature", "Heatwave", "CDHW"),
 col = c("grey60", "red", "darkred"),
 pch = c(NA, 16, 16),
 lty = c(1, NA, NA),
 pt.cex = 0.7,
 bty = "n",
 cex = 0.8,
 xpd = TRUE
)

dev.off()


## 2022
tiff(
 filename = "CDHW_2022.tiff",
 width = 2000,
 height = 1200,
 res = 330
)

plot(
 tmax_2022$date,
 tmax_2022$Tmax,
 type = "l",
 col = "grey60",
 xlab = "",
 ylab = "Daily maximum temperature (°C)",
 bty  = "l",
 xaxt = "n"
)

# Custom x-axis for months
axis.Date(
 side = 1,
 at = seq(as.Date("2022-01-01"), as.Date("2022-12-01"), by = "month"),
 format = "%b"
)
mtext("2022", side = 1, line = 2.5)

# Heatwave days
points(
 tmax_2022$date[tmax_2022$heatwave_event_day],
 tmax_2022$Tmax[tmax_2022$heatwave_event_day],
 col = "red",
 pch = 16,
 cex = 0.5
)

# CDHW days
points(
 tmax_2022$date[tmax_2022$CDHW_event_day],
 tmax_2022$Tmax[tmax_2022$CDHW_event_day],
 col = "darkred",
 pch = 16,
 cex = 0.5
)

# Legend
legend(
 "topright",
 legend = c("Daily max temperature", "Heatwave", "CDHW"),
 col = c("grey60", "red", "darkred"),
 pch = c(NA, 16, 16),
 lty = c(1, NA, NA),
 pt.cex = 0.7,
 bty = "n",
 cex = 0.8
)

dev.off()


##### 
# 2022 Dual-axis plot with SPEI-3 threshold
par(mar = c(5, 5, 4, 5))  # bottom, left, top, right

# Left axis: daily Tmax
plot(
 tmax_2022$date,
 tmax_2022$Tmax,
 type = "l",
 col = "grey60",
 xlab = "",
 ylab = "Daily maximum temperature (°C)",
 bty = "l",
 xaxt = "n",
 lwd = 1.5,
 cex.lab = 1.5,    # bigger axis labels
 cex.axis = 1.3    # bigger tick labels
)

# Custom x-axis for months
axis.Date(
 side = 1,
 at = seq(as.Date("2022-01-01"), as.Date("2022-12-01"), by = "month"),
 format = "%b",
 cex.axis = 1.3
)
mtext("2022", side = 1, line = 2.5, cex = 1.5)

# Heatwave days
points(
 tmax_2022$date[tmax_2022$heatwave_event_day],
 tmax_2022$Tmax[tmax_2022$heatwave_event_day],
 col = "red",
 pch = 16,
 cex = 0.5
)

# CDHW days
points(
 tmax_2022$date[tmax_2022$CDHW_event_day],
 tmax_2022$Tmax[tmax_2022$CDHW_event_day],
 col = "darkred",
 pch = 16,
 cex = 0.5
)

# Adding right axis: SPEI-3
par(new = TRUE)  # overlay second plot
plot(
 tmax_2022$date,
 tmax_2022$SPEI3,
 type = "l",
 col = "#0047AB",
 axes = FALSE,
 xlab = "",
 ylab = "",
 lwd = 2,
 ylim = c(min(tmax_2022$SPEI3, -2), max(tmax_2022$SPEI3, 2))
)
axis(side = 4, col.axis = "black", col = "black", cex.axis = 1.3)
mtext("SPEI-3", side = 4, line = 2, col = "black", cex = 1.5)

# Horizontal line for drought threshold
abline(h = -1, col = "black", lty = 2, lwd = 1)

# Legend
legend(
 "topright",
 legend = c("Daily max temperature", "Heatwave", "CDHW", "SPEI-3", "Drought threshold (SPEI-3 < -1)"),
 col = c("grey60", "red", "darkred", "darkblue", "darkblue"),
 pch = c(NA, 16, 16, NA, NA),
 lty = c(1, NA, NA, 1, 2),
 pt.cex = 0.7,
 bty = "n",
 cex = 0.8
)


######## CDHW intensity–duration scatterplots

# Prepare CDHW metrics
CDHW_plot_df <- CDHW_events %>%
 mutate(
  intensity = severity_heat / duration,
  year = year(start_date)
 )

# Scatterplot: Duration vs Intensity
ggplot(CDHW_plot_df,
       aes(x = duration, y = intensity, color = year)) +
 geom_point(size = 3) +
 scale_color_viridis_c() +
 geom_smooth(method = "lm", se = FALSE, col = "grey40") +
 labs(
  title = "CDHW Intensity–Duration Relationship (Sussex)",
  x = "Duration (days)",
  y = "Mean heatwave intensity (°C above threshold)",
  color = "Year"
 ) +
 theme_minimal(base_size = 14)


############## CDHW Severity Classification
tmax_df <- tmax_df %>%
 mutate(
  Tmax_std = scale(Tmax)  # z-score: (Tmax - mean)/sd
 )

tmax_df <- tmax_df %>%
 mutate(
  daily_CDHW_severity = Tmax_std * (-SPEI3)  # higher values = more severe compound
 )

CDHW_events <- tmax_df %>%
 filter(!is.na(hw_event_id) & CDHW_day) %>%
 group_by(hw_event_id) %>%
 summarise(
  start_date = min(date),
  end_date   = max(date),
  duration   = n(),
  Tmax_mean  = mean(Tmax),
  Tmax_max   = max(Tmax),
  mean_SPEI3 = mean(SPEI3),
  event_CDHW_severity = sum(daily_CDHW_severity),
  .groups = "drop"
 )

CDHW_events <- CDHW_events %>%
 mutate(
  severity_norm = scale(event_CDHW_severity),
  CDHW_class = cut(
   severity_norm,
   breaks = quantile(severity_norm, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
   labels = c("Low", "Moderate", "Severe", "Extreme"),
   include.lowest = TRUE
  )
 )

write.csv(
 CDHW_events,
 "CDHW_events.csv",
 row.names = FALSE
)

################ CDHW in 2022
# Prepare CDHW_JJA_2022 with a "type" column for legend
CDHW_JJA_2022 <- tmax_df %>%
 filter(month(date) %in% 6:8, year(date) == 2022) %>%
 left_join(
  CDHW_events %>% select(hw_event_id, CDHW_class),
  by = c("hw_event_id")
 ) %>%
 mutate(
  point_type = case_when(
   CDHW_event_day ~ as.character(CDHW_class),        # CDHW colored by severity
   heatwave_event_day ~ "Heatwave",                  # Heatwave-only
   TRUE ~ NA_character_
  ),
  # Factor for legend order
  point_type = factor(point_type, levels = c("Heatwave", "Low", "Moderate", "Severe", "Extreme"))
 )

# Define colors
point_colors <- c(
 "Low"      = "#E6C200",
 "Moderate" = "orange",
 "Severe"   = "#CC0000",
 "Extreme"  = "#660000"
)

# 15-day breaks for x-axis
day_breaks <- seq(as.Date("2022-06-01"), as.Date("2022-08-31"), by = "15 days")

# Plot
p <- ggplot(CDHW_JJA_2022, aes(x = date)) +
 geom_line(aes(y = Tmax), color = "grey60") +
 
 # Points (heatwave + CDHW)
 geom_point(
  data = CDHW_JJA_2022 %>% filter(!is.na(point_type)),
  aes(y = Tmax, color = point_type),
  size = 2
 ) +
 
 scale_color_manual(
  name = "Event type / Severity",
  values = point_colors
 ) +
 
 scale_x_date(
  breaks = day_breaks,
  labels = function(x) format(x, "%d %B") 
 ) +
 
 labs(
  x = "2022",
  y = "Daily maximum temperature (°C)"
 ) +
 
 theme_minimal(base_size = 14) +
 theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.5),
  axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 14),
  axis.text.y = element_text(angle = 0, hjust = 0.5, color = "black", size = 14),
  legend.position = "right"
 )

# Print plot
p

# Save plot as TIFF
ggsave(
 filename = "CDHW_JJA_2022.tif",
 plot = p,
 device = "tiff",
 dpi = 330,
 width = 25,  
 height = 12,
 units = "cm"
)


# Prepare CDHW_JJA_2022 with severity classification
CDHW_JJA_2022 <- tmax_df %>%
 filter(month(date) %in% 6:8, year(date) == 2022) %>%
 left_join(
  CDHW_events %>% select(hw_event_id, CDHW_class),
  by = c("hw_event_id")
 ) %>%
 mutate(
  point_type = case_when(
   CDHW_event_day ~ as.character(CDHW_class),
   heatwave_event_day ~ "Heatwave",
   TRUE ~ NA_character_
  ),
  point_type = factor(
   point_type,
   levels = c("Heatwave", "Low", "Moderate", "Severe", "Extreme")
  )
 )

# Define colors
point_colors <- c(
 "Heatwave" = "red",
 "Low"      = "#E6C200",
 "Moderate" = "orange",
 "Severe"   = "#CC0000",
 "Extreme"  = "#660000"
)

# Define axis breaks

day_breaks <- seq(as.Date("2022-06-01"),
                  as.Date("2022-08-31"),
                  by = "15 days")

# Rescale SPEI to match Tmax axis (required for ggplot dual axis)

tmax_min <- min(CDHW_JJA_2022$Tmax, na.rm = TRUE)
tmax_max <- max(CDHW_JJA_2022$Tmax, na.rm = TRUE)

spei_min <- -2.5
spei_max <-  2.5

scale_spei_to_tmax <- function(x) {
 (x - spei_min) / (spei_max - spei_min) *
  (tmax_max - tmax_min) + tmax_min
}

scale_tmax_to_spei <- function(x) {
 (x - tmax_min) / (tmax_max - tmax_min) *
  (spei_max - spei_min) + spei_min
}

# Plot
p <- ggplot(CDHW_JJA_2022, aes(x = date)) +
 
 # Tmax line
 geom_line(aes(y = Tmax),
           color = "grey60",
           linewidth = 1) +
 
 # SPEI-3 line (rescaled)
 geom_line(aes(y = scale_spei_to_tmax(SPEI3)),
           color = "#0047AB",
           linewidth = 1.3) +
 
 # Drought threshold (SPEI = -1)
 geom_hline(
  yintercept = scale_spei_to_tmax(-1),
  linetype = "dashed",
  color = "black",
  linewidth = 0.8
 ) +
 
 # Event points
 geom_point(
  data = CDHW_JJA_2022 %>% filter(!is.na(point_type)),
  aes(y = Tmax, color = point_type),
  size = 2
 ) +
 
 scale_color_manual(
  name = "Event type / Severity",
  values = point_colors
 ) +
 
 scale_x_date(
  breaks = day_breaks,
  labels = function(x) format(x, "%d %B")
 ) +
 
 scale_y_continuous(
  name = "Daily maximum temperature (°C)",
  sec.axis = sec_axis(
   trans = ~ scale_tmax_to_spei(.),
   name = "SPEI-3"
  )
 ) +
 
 labs(x = "2022") +
 
 theme_minimal(base_size = 14) +
 theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.5),
  axis.text.x = element_text(color = "black", size = 14),
  axis.text.y = element_text(color = "black", size = 14),
  axis.title.y.right = element_text(color = "#0047AB"),
  legend.position = "right"
 )

# Print plot
p
