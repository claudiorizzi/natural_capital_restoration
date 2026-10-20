# -*- coding: utf-8 -*-
"""
04_maps.py - Hexbin map of TREES farm locations (Figure: farm_location_map).

Changes vs v1
-------------
- Input switched to the EVI estimation panel `final/EVI_farms_50_km_longer_time.dta`
  (v1 used `final/EVI_farms_longer_time.dta`).
- The sample is cleaned exactly as at the top of `evi_analysis` in
  `03_main_estimates.do`:
      gen event_window = year - calyear_treat if calyear_treat != 0
      keep if (event_window >= -3 & event_window <= 4) | calyear_treat == 0
- Only treated farms (calyear_treat != 0) are mapped, matching the caption
  ("farms that participated in the TREES program"). v1 filtered on IDs containing
  "FG", but 4,166 FG farms in the new panel are never treated (calyear_treat == 0).
- New: each farm is assigned to the Natural Earth country polygon it falls in and
  the farms that the map places in Ethiopia are exported to `Data/Processed/`.

Inputs  (DATA_DIR from .env)
    final/EVI_farms_pseudo.dta   (the paper's file is EVI_farms_50_km_longer_time.dta)
    country_boundaries/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp
Outputs
    OVERLEAF_FIGURES/farm_location_map.png
    Data/Processed/farms_mapped_in_ethiopia.csv  (+ .dta)
"""

# Replication package: identical to the authors' Code/6_1_maps_v2.py except that the
# paths read from .env are resolved relative to the package root (so the
# shipped .env with relative paths works from any working directory) and the Ethiopia check file goes to output/processed.

import os
from pathlib import Path

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")  # non-interactive backend: safe under run_all.py

# =======================================================
# Configuration
# =======================================================
# Load paths from project-root .env so they are configurable per machine
PROJECT_ROOT = Path(__file__).resolve().parents[1]
load_dotenv(PROJECT_ROOT / ".env")
DATA_DIR = (PROJECT_ROOT / os.environ["DATA_DIR"]).resolve()
FARMS_PATH = DATA_DIR / "final/EVI_farms_pseudo.dta"
MAP_PATH = DATA_DIR / "country_boundaries/ne_110m_admin_0_countries"
OUTPUT_DIR = (PROJECT_ROOT / os.environ["OVERLEAF_FIGURES"]).resolve()
PROCESSED_DIR = PROJECT_ROOT / "output" / "processed"

# Only the columns the map needs: the .dta is ~230 MB / 118 columns
COLS = [
    "forestgardenid", "year", "calyear_treat", "fe_cell", "lat", "lon",
    "iso_trees", "projectunit", "proj_area", "status",
]


# =======================================================
# 1. Load and clean the panel as in evi_analysis (03_main_estimates.do)
# =======================================================
print(f"Loading {FARMS_PATH} ...")
df = pd.read_stata(FARMS_PATH, columns=COLS, convert_categoricals=False)
print(f"  {len(df):,} farm-years, {df['forestgardenid'].nunique():,} farms")

# gen event_window = year - calyear_treat if calyear_treat != 0
df["event_window"] = np.where(df["calyear_treat"] != 0,
                              df["year"] - df["calyear_treat"], np.nan)

# keep if (event_window >= -3 & event_window <= 4) | calyear_treat == 0
# (NaN comparisons are False in pandas, as the Stata condition drops missing
#  calyear_treat rows too, since missing <= 4 is false in Stata.)
in_window = df["event_window"].between(-3, 4)
df = df[in_window | (df["calyear_treat"] == 0)].copy()

# egen fe_cellyr = group(fe_cell year)  -- not used by the map, kept for parity
df["fe_cellyr"] = df.groupby(["fe_cell", "year"]).ngroup() + 1

print(f"  after window: {len(df):,} farm-years, {df['forestgardenid'].nunique():,} farms")

# =======================================================
# 2. Farm-level sample for the map: treated farms with coordinates
# =======================================================
farms = (df.sort_values(["forestgardenid", "year"])
           .drop_duplicates(subset="forestgardenid")
           .copy())
n_fg_never = farms.loc[farms["forestgardenid"].str.contains("FG", na=False)
                       & (farms["calyear_treat"] == 0)].shape[0]
farms = farms[farms["calyear_treat"] != 0].copy()          # treated only
n_nocoord = farms["lat"].isna().sum() + farms["lon"].isna().sum()
farms = farms.dropna(subset=["lat", "lon"]).copy()
print(f"  treated farms mapped: {len(farms):,} "
      f"(never-treated FG farms excluded: {n_fg_never:,}; "
      f"treated farms without coordinates: {n_nocoord})")

gdf = gpd.GeoDataFrame(
    farms, geometry=gpd.points_from_xy(farms["lon"], farms["lat"]), crs="EPSG:4326"
)

# =======================================================
# 3. Natural Earth boundaries + country the map assigns to each farm
# =======================================================
africa = gpd.read_file(MAP_PATH / "ne_110m_admin_0_countries.shp")
africa = africa[africa["CONTINENT"] == "Africa"]

# Point-in-polygon on the same (coarse, 1:110m) polygons the map draws, so this
# is literally "the country the map thinks the farm is in".
joined = gpd.sjoin(
    gdf, africa[["NAME", "ADM0_A3", "geometry"]], how="left", predicate="within"
).rename(columns={"NAME": "map_country", "ADM0_A3": "map_iso3"})
joined = joined.drop(columns="index_right")
joined["map_iso3"] = joined["map_iso3"].fillna("NONE")   # off every polygon (sea / gaps)

mismatch = joined[joined["map_iso3"] != joined["iso_trees"]]
print("\nFarms whose map country differs from iso_trees (TREES country):")
print(mismatch.groupby(["iso_trees", "map_iso3"]).size().to_string())

# Export farms the map places in Ethiopia (TREES does not operate there)
eth = joined[joined["map_iso3"] == "ETH"][
    ["forestgardenid", "lat", "lon", "iso_trees", "projectunit", "proj_area",
     "status", "calyear_treat", "map_country", "map_iso3"]
].sort_values("forestgardenid")
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
eth.to_csv(PROCESSED_DIR / "farms_mapped_in_ethiopia.csv", index=False)
eth.to_stata(PROCESSED_DIR / "farms_mapped_in_ethiopia.dta",
             write_index=False, version=118)
print(f"\n{len(eth)} farms mapped in Ethiopia -> "
      f"{PROCESSED_DIR / 'farms_mapped_in_ethiopia.csv'} (+ .dta)")
if len(eth):
    print(eth[["iso_trees", "projectunit"]].value_counts().to_string())

# =======================================================
# 4. Plot hexbin map
# =======================================================
# Exclude the farms the map places in Ethiopia: 5 Tanzania CDI 1 farms whose
# latitude sign is flipped (7.58 instead of -7.58). Decision: Claudio, 2026-09-25.
gdf = gdf[~gdf["forestgardenid"].isin(eth["forestgardenid"])].copy()
# Remove obvious bad coordinates outside Africa's bounding box
gdf = gdf[gdf["lat"].between(-40, 40) & gdf["lon"].between(-25, 60)].copy()
print(f"Farms plotted after excluding Ethiopia-mapped farms: {len(gdf):,}")

fig, ax = plt.subplots(figsize=(7, 7))

# Plot Africa boundary
africa.boundary.plot(ax=ax, linewidth=0.6, color="black")

# Medium-light green -> dark green
cmap = LinearSegmentedColormap.from_list("green_gradient", ["#66c2a4", "#006d2c"])

hb = ax.hexbin(gdf.geometry.x, gdf.geometry.y, gridsize=40, mincnt=1,
               cmap=cmap, linewidths=0)

# Colorbar
cb = fig.colorbar(hb, ax=ax, shrink=0.5, aspect=20, pad=0.02)
cb.set_label("Number of farms per grid")

# Formatting
ax.set_title("")
ax.set_xlim(-20, 55)
ax.set_ylim(-35, 20)
ax.axis("off")
plt.tight_layout()

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
out_png = OUTPUT_DIR / "farm_location_map.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out_png}")
