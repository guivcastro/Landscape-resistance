"""
NDVI resistance using ARIMA per habitat type

Baseline: 1991–2020
Expected: NDVI_2022 using data through 2021

Filenames format:
NDVI_broadleaf_1991.tif
"""
import statsmodels
import numpy as np
import rasterio
from glob import glob
import os
from collections import defaultdict
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# Data directory
DATA_DIR = '/work/scratch-pw5/guicas/NDVI/Resistance/raw_habitat'
OUTPUT_DIR = "/work/scratch-pw5/guicas/NDVI/Resistance/outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASELINE_START = 1991
BASELINE_END   = 2020
LAG_YEAR       = 2021
TARGET_YEAR    = 2022

Z_VALUE = 1.96

# Find files
files = sorted(glob(os.path.join(DATA_DIR, "*.tif")))

if len(files) == 0:
    raise FileNotFoundError("No rasters found.")

print(f"Found {len(files)} rasters")

# Group files by habitat
habitat_files = defaultdict(list)

for f in files:

    name = os.path.basename(f).replace(".tif", "")
    name = name.replace("NDVI_", "")

    habitat, year = name.rsplit("_", 1)
    habitat_files[habitat].append((int(year), f))

print("\nHabitats detected:")
for h in habitat_files:
    print(" •", h)

# Prepare output arrays
meta = None
combined_expected = None
combined_resistance = None
combined_significance = None
combined_process_sd = None

def write_tif(filename, array, meta, dtype="float32"):
    path = os.path.join(OUTPUT_DIR, filename)

    m = meta.copy()
    m.update(dtype=dtype, count=1, compress="lzw")

    with rasterio.open(path, "w", **m) as dst:
        dst.write(array.astype(dtype), 1)

summary_rows = []

# Loop through habitat types
for habitat, year_files in habitat_files.items():

    print(f"\nProcessing habitat: {habitat}")

    # Sort files by year
    year_files = sorted(year_files)
    years = np.array([yf[0] for yf in year_files])

    rasters = []

    for _, f in year_files:
        with rasterio.open(f) as src:
            rasters.append(src.read(1))

            if meta is None:
                meta = src.meta.copy()

    ndvi = np.stack(rasters)  # (time, rows, cols)

    # Initialise combined rasters once
    if combined_expected is None:
        shape = ndvi.shape[1:]

        combined_expected = np.full(shape, np.nan, dtype="float32")
        combined_resistance = np.full(shape, np.nan, dtype="float32")
        combined_significance = np.full(shape, -2, dtype="int8")
        combined_process_sd = np.full(shape, np.nan, dtype="float32")

    # Habitat-specific output rasters
    shape = ndvi.shape[1:]

    hab_expected = np.full(shape, np.nan, dtype="float32")
    hab_resistance = np.full(shape, np.nan, dtype="float32")
    hab_significance = np.full(shape, -2, dtype="int8")
    hab_process_sd = np.full(shape, np.nan, dtype="float32")

    # Build habitat timeseries
    needed_years = (years >= BASELINE_START) & (years <= LAG_YEAR)

    if np.sum(needed_years) < 10:
        print("Too few years — skipping habitat")
        continue

    habitat_series = np.nanmean(ndvi[needed_years], axis=(1, 2))

    # remove any NaNs in the series
    valid_idx = np.isfinite(habitat_series)
    habitat_series = habitat_series[valid_idx]

    if len(habitat_series) < 10:
        print("Insufficient valid data — skipping habitat")
        continue

    # Stationarity check using ADF test
    adf_result = adfuller(habitat_series)
    adf_stat = adf_result[0]
    adf_pval = adf_result[1]

    if adf_pval < 0.05:
        print(f"ADF test for {habitat}: stationary (p={adf_pval:.3f})")
        d_candidates = [0]
    else:
        print(f"ADF test for {habitat}: NOT stationary (p={adf_pval:.3f})")
        d_candidates = [1]

    # Auto ARIMA search using AICc + parsimony
    models = []

    for p in range(3):
        for d in d_candidates:
            for q in range(3):

                try:
                    result = ARIMA(habitat_series, order=(p, d, q)).fit()

                    aic = result.aic
                    k = result.params.size
                    n = result.nobs

                    # Avoid division errors
                    if (n - k - 1) <= 0:
                        continue

                    aicc = aic + (2 * k * (k + 1)) / (n - k - 1)

                    models.append({
                        "aicc": aicc,
                        "order": (p, d, q),
                        "result": result,
                        "complexity": p + q   # measure of parsimony
                    })

                except:
                    continue

    if len(models) == 0:
        print("ARIMA failed — skipping habitat")
        continue

    # Sort by AICc
    models = sorted(models, key=lambda x: x["aicc"])

    # Print AICc table for this habitat
    print(f"\nAICc values for habitat: {habitat}")
    print(f"{'Model':>10} | {'AICc':>10}")
    print("-"*23)
    for m in models:
        print(f"ARIMA{m['order']} | {m['aicc']:10.2f}")

    # Select the single best model (lowest AICc)
    best = models[0]

    best_model = best["result"]
    best_order = best["order"]

    print(f"Selected best model: ARIMA{best_order} (AICc={best['aicc']:.2f})")

    from statsmodels.stats.diagnostic import acorr_ljungbox

    # Check residual autocorrelation (white noise)
    lb = acorr_ljungbox(best_model.resid, lags=[10], return_df=True)
    lb_stat = lb['lb_stat'].values[0]
    lb_pval = lb['lb_pvalue'].values[0]

    if lb_pval > 0.05:
        print(f"Ljung-Box test: residuals are white noise (p={lb_pval:.3f})")
    else:
        print(f"Ljung-Box test WARNING: residuals may have autocorrelation (p={lb_pval:.3f})")

    from statsmodels.stats.stattools import jarque_bera

    # Jarque-Bera test for residual normality
    jb_stat, jb_p, _, _ = jarque_bera(best_model.resid)
    if jb_p > 0.05:
        print(f"Jarque-Bera test: residuals are approximately normal (p={jb_p:.3f})")
    else:
        print(f"Jarque-Bera WARNING: residuals deviate from normality (p={jb_p:.3f})")
  
    # Forecast 2022
    forecast = best_model.get_forecast(steps=1)
    expected_mean = forecast.predicted_mean[0]
    se_2022 = forecast.se_mean[0]
    process_sd = se_2022

    # Load 2021 and 2022 raster data
    try:
        ndvi_2021 = ndvi[years == LAG_YEAR][0]
        observed_2022 = ndvi[years == TARGET_YEAR][0]
    except:
        print("Missing 2021 or 2022 — skipping habitat")
        continue

    # Compute expected 2022 from 2021 pixels
    mean_2021 = np.nanmean(ndvi_2021)
    if mean_2021 < 1e-6:
        print("Habitat mean near zero — skipping")
        continue

    scaling_ratio = expected_mean / mean_2021
    expected_2022 = ndvi_2021 * scaling_ratio
    expected_2022[expected_2022 < 0] = np.nan  # remove negative expected NDVI

    # Compute resistance and filter
    resistance = observed_2022 / expected_2022

    # Remove negative resistance
    resistance[resistance < 0] = np.nan

    # Remove extreme resistance values 98th percentile
    lower_res = np.nanpercentile(resistance, 1)
    upper_res = np.nanpercentile(resistance, 99)
    resistance[(resistance < lower_res) | (resistance > upper_res)] = np.nan

    # Confidence interval
    lower_ci = expected_2022 - Z_VALUE * se_2022
    upper_ci = expected_2022 + Z_VALUE * se_2022

    significance = np.zeros_like(observed_2022, dtype=np.int8)
    significance[observed_2022 < lower_ci] = -1
    significance[observed_2022 > upper_ci] = 1

    # Merge into combined rasters
    habitat_mask = np.isfinite(ndvi_2021)
    combined_expected[habitat_mask] = expected_2022[habitat_mask]
    combined_resistance[habitat_mask] = resistance[habitat_mask]
    combined_significance[habitat_mask] = significance[habitat_mask]
    combined_process_sd[habitat_mask] = process_sd

    # Habitats
    hab_expected[habitat_mask] = expected_2022[habitat_mask]
    hab_resistance[habitat_mask] = resistance[habitat_mask]
    hab_significance[habitat_mask] = significance[habitat_mask]
    hab_process_sd[habitat_mask] = process_sd

    valid_mask = hab_significance != -2   # ignore untouched pixels

    total_pixels = np.sum(valid_mask)

    positive_pixels = np.sum(hab_significance == 1)
    negative_pixels = np.sum(hab_significance == -1)
    nonsig_pixels  = np.sum(hab_significance == 0)
  
    if total_pixels > 0:
        pos_pct = (positive_pixels / total_pixels) * 100
        neg_pct = (negative_pixels / total_pixels) * 100
        nonsig_pct = (nonsig_pixels / total_pixels) * 100
    else:
        pos_pct = neg_pct = nonsig_pct = np.nan

    summary_rows.append({
        "Habitat": habitat,
        "Total_pixels": total_pixels,
        "Positive_pixels": positive_pixels,
        "Negative_pixels": negative_pixels,
        "Non_significant_pixels": nonsig_pixels,
        "Positive_%": pos_pct,
        "Negative_%": neg_pct,
        "Non_significant_%": nonsig_pct,
        "Ljung_Box_p_lag10": lb_pval
    })

    # EXPORT HABITAT HERE
    safe_habitat = habitat.replace(" ", "_")

    hab_ci_range = Z_VALUE * hab_process_sd * 2

    sig_mask_h = hab_significance != 0
    pos_mask_h = hab_significance == 1
    neg_mask_h = hab_significance == -1
    nonsig_mask_h = hab_significance == 0

    write_tif(f"NDVI_expected_2022_{safe_habitat}.tif", hab_expected, meta)
    write_tif(f"NDVI_resistance_2022_{safe_habitat}.tif", hab_resistance, meta)
    write_tif(f"NDVI_significance_2022_{safe_habitat}.tif", hab_significance, meta, dtype="int8")
    write_tif(f"NDVI_process_sd_{safe_habitat}.tif", hab_process_sd, meta)

    write_tif(f"NDVI_CI_range_2022_{safe_habitat}.tif", hab_ci_range, meta)

    write_tif(
        f"Resistance_significant_2022_{safe_habitat}.tif",
        np.where(sig_mask_h, hab_resistance, np.nan),
        meta
    )

    write_tif(
        f"Resistance_positive_2022_{safe_habitat}.tif",
        np.where(pos_mask_h, hab_resistance, np.nan),
        meta
    )

    write_tif(
        f"Resistance_negative_2022_{safe_habitat}.tif",
        np.where(neg_mask_h, hab_resistance, np.nan),
        meta
    )

    write_tif(
        f"Resistance_nonsignificant_2022_{safe_habitat}.tif",
        np.where(nonsig_mask_h, hab_resistance, np.nan),
        meta
    )

# Export outputs
write_tif("NDVI_expected_2022_combined.tif", combined_expected, meta)
write_tif("NDVI_resistance_2022_combined.tif", combined_resistance, meta)
write_tif("NDVI_significance_2022_combined.tif", combined_significance, meta, dtype="int8")
write_tif("NDVI_process_sd_combined.tif", combined_process_sd, meta)

# CI width (prediction uncertainty)
combined_ci_range = Z_VALUE * combined_process_sd * 2

sig_mask = combined_significance != 0
pos_mask = combined_significance == 1
neg_mask = combined_significance == -1
nonsig_mask = combined_significance == 0

write_tif(
    "NDVI_CI_range_2022_combined.tif",
    combined_ci_range,
    meta
)

write_tif(
    "Resistance_significant_2022_combined.tif",
    np.where(sig_mask, combined_resistance, np.nan),
    meta
)

write_tif(
    "Resistance_positive_2022_combined.tif",
    np.where(pos_mask, combined_resistance, np.nan),
    meta
)

write_tif(
    "Resistance_negative_2022_combined.tif",
    np.where(neg_mask, combined_resistance, np.nan),
    meta
)

write_tif(
    "Resistance_nonsignificant_2022_combined.tif",
    np.where(nonsig_mask, combined_resistance, np.nan),
    meta
)

print("\nHabitat ARIMA resistance rasters created")

valid_mask = combined_significance != -2

total_pixels = np.sum(valid_mask)

positive_pixels = np.sum(combined_significance == 1)
negative_pixels = np.sum(combined_significance == -1)
nonsig_pixels  = np.sum(combined_significance == 0)

summary_rows.append({
    "Habitat": "COMBINED",
    "Total_pixels": total_pixels,
    "Positive_pixels": positive_pixels,
    "Negative_pixels": negative_pixels,
    "Non_significant_pixels": nonsig_pixels,
    "Positive_%": (positive_pixels / total_pixels) * 100,
    "Negative_%": (negative_pixels / total_pixels) * 100,
    "Non_significant_%": (nonsig_pixels / total_pixels) * 100
})

summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    os.path.join(OUTPUT_DIR, "NDVI_resistance_pixel_summary.csv"),
    index=False,
    float_format="%.2f"
)

print("\nPixel summary exported → NDVI_resistance_pixel_summary.csv")
