# -*- coding: utf-8 -*-
"""
01_make_pseudo_data.py - generate the SYNTHETIC demo datasets shipped with this
replication package.

Why synthetic data
------------------
The household-survey and farm-location data used in the paper were shared by
Trees for the Future (TREES) under a data-use agreement. They describe
individual rural households (household composition, food security, livestock,
farm size) and carry farm-level GPS coordinates, so the records cannot be
de-identified without destroying the information the analysis relies on. The
raw and processed data therefore cannot be redistributed. This script creates
pseudo datasets with the SAME variable names, types, panel structure and
sample-definition flags as the processed files read by the Stata scripts, so
that every program in `02_summary_stats.do` and `03_main_estimates.do`
(and the map script `04_maps.py`) runs end-to-end.

The two public inputs of the paper (the World Bank LSMS-ISA Tanzania National
Panel Survey 2020-21 and the local-GDP grid of Rossi-Hansberg and Zhang; see
the README for the download pages) are not shipped either: for simplicity of
the replication code the package contains no raw files at all, and the pseudo
datasets carry synthetic versions of the variables derived from them.

What the pseudo data are, and are not
-------------------------------------
* Nothing in these files is derived from a real farm. Farm IDs, coordinates,
  project names and every outcome are drawn from random-number generators; the
  moments (means, dispersion, shares, cohort sizes) are only loosely matched to
  the real data so the demo produces sensible-looking output. No identifier,
  name or coordinate of a real farm or household appears in the data or
  anywhere in the code.
* The pseudo data are a synthetic version of the data, NOT a copy: the numbers
  of farms, surveys and farm-years (1,600 farms / 2,000 surveys; 5,000 farms x
  8 years = 40,000 farm-years) are set for a quick demo and do not match the
  samples of the paper.
  Point estimates obtained on the pseudo data are NOT the paper's estimates.
* The household panel is written with ONLY the variables used by the paper's
  scripts (the list KEEP_HH at the end of section 1; 88 variables). The
  data-generating process below simulates the raw-survey structure in more
  detail (intermediate string answers, unused characteristics, extra
  fixed-effect ids) because the kept variables are built from it, but these
  intermediate variables are not written to the file.
* Treatment effects built into the data-generating process: after the first
  survey, livestock ownership +6 percentage points and herd size +25 log
  points; trees +50 log points and crops +15 log points; small improvements in
  food security; and a dynamic EVI effect of 0.3, 1.0, 1.5, 2.0 and 2.5 EVI
  points (x 100) in years 0 to 4 after treatment, with no pre-trend.

Files produced (all under ../data/)
-----------------------------------
final/TREES_data_pseudo.dta   household panel, restricted to the variables the scripts use (KEEP_HH)
final/EVI_farms_pseudo.dta    satellite (EVI) farm-year panel
final/tanzania_merged.dta     stand-in for the World Bank LSMS-ISA Tanzania benchmark

Usage (from the repository root, with uv):
    uv run python code/01_make_pseudo_data.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261019
rng = np.random.default_rng(SEED)

ROOT = Path(__file__).resolve().parents[1]          # repository root
DATA = ROOT / "data"
FINAL = DATA / "final"
FINAL.mkdir(parents=True, exist_ok=True)

# =======================================================
# Shared geography: countries, fictional project areas, GDP grid cells
# =======================================================
# ISO code, country name, share of household farms, share of satellite farms,
# bounding box used to scatter fictional farm locations (lat_min, lat_max,
# lon_min, lon_max). The boxes sit inside each country's interior so the map
# demo places points on land; they are NOT real project locations.
COUNTRIES = {
    "KEN": dict(name="Kenya",    hh=0.25, sat=0.36, box=(-1.5, 1.0, 34.3, 35.8)),
    "SEN": dict(name="Senegal",  hh=0.30, sat=0.30, box=(14.2, 15.3, -16.4, -14.5)),
    "TZA": dict(name="Tanzania", hh=0.23, sat=0.15, box=(-7.0, -4.0, 33.5, 36.0)),
    "UGA": dict(name="Uganda",   hh=0.22, sat=0.12, box=(0.3, 1.5, 32.8, 34.0)),
    "MLI": dict(name="Mali",     hh=0.00, sat=0.07, box=(12.0, 13.5, -8.5, -6.0)),
}
N_AREAS = {"KEN": 7, "SEN": 7, "TZA": 5, "UGA": 5, "MLI": 3}   # project areas per country
UNITS_PER_AREA = 4                                              # project units per area

# project areas / units / grid cells ------------------------------------
areas = []
for iso, k in N_AREAS.items():
    for a in range(1, k + 1):
        lat0 = rng.uniform(COUNTRIES[iso]["box"][0], COUNTRIES[iso]["box"][1])
        lon0 = rng.uniform(COUNTRIES[iso]["box"][2], COUNTRIES[iso]["box"][3])
        areas.append(dict(iso=iso, proj_area=f"{COUNTRIES[iso]['name']}: Area {a:02d}",
                          lat0=lat0, lon0=lon0,
                          # each area spans 1-2 GDP grid cells (0.25 degree)
                          cells=[int(27000 + 1000 * list(N_AREAS).index(iso) + 20 * a + j)
                                 for j in range(rng.integers(1, 3))]))
AREAS = pd.DataFrame(areas)
UNITS = pd.DataFrame([dict(iso=r.iso, proj_area=r.proj_area, projectunit=f"{r.proj_area} {u}",
                           lat0=r.lat0, lon0=r.lon0, cell_id=int(rng.choice(r.cells)),
                           # project-unit start year: when the TREES project in that unit began
                           unit_start=int(rng.choice([2014, 2015, 2016, 2018, 2019, 2020, 2021, 2022, 2023, 2024],
                                                     p=[.02, .02, .03, .08, .09, .16, .12, .05, .19, .24])))
                      for r in AREAS.itertuples() for u in range(1, UNITS_PER_AREA + 1)])
ALL_CELLS = sorted(set(UNITS.cell_id))

# cell x year local GDP (Rossi-Hansberg & Zhang style): gcp in bn 2017 USD,
# population per cell, GDP per capita. Mild common growth + cell-specific noise.
YEARS = list(range(2014, 2025))
cell_rows = []
for c in ALL_CELLS:
    base_gcp = np.exp(rng.normal(np.log(0.15), 1.0))           # heavy right tail like the real grid
    pop = np.exp(rng.normal(np.log(150_000), 0.9))
    g = rng.normal(0.02, 0.01)
    for t, y in enumerate(YEARS):
        gcp = base_gcp * np.exp(g * t + rng.normal(0, 0.04))
        cell_rows.append(dict(cell_id=c, year=y, gcp_2017usd=gcp, pop_cell=pop * (1 + 0.02 * t),
                              gcpc_2017usd=gcp / (pop * (1 + 0.02 * t))))
CELLS = pd.DataFrame(cell_rows).sort_values(["cell_id", "year"])
CELLS["l1_gcpc_2017usd"] = CELLS.groupby("cell_id").gcpc_2017usd.shift(1)
CELLS["ch_gcpc_2017usd"] = (CELLS.gcpc_2017usd - CELLS.l1_gcpc_2017usd) / CELLS.l1_gcpc_2017usd


def winsor(x, lo=1, hi=99):
    """winsor2-style winsorisation at the lo/hi percentiles (NaN kept)."""
    x = pd.Series(x, dtype="float64")
    if x.notna().sum() == 0:
        return x
    a, b = np.nanpercentile(x, [lo, hi])
    return x.clip(lower=a, upper=b)


def group_id(df, cols):
    """Stata egen group(): consecutive integer ids by sorted unique combination."""
    return df.groupby(cols, sort=True, dropna=False).ngroup() + 1


def yes_no(p, n, missing=0.0):
    s = np.where(rng.random(n) < p, "Yes", "No").astype(object)
    s[rng.random(n) < missing] = ""
    return s


# =======================================================
# 1. Household panel  (final/TREES_data_pseudo.dta)
# =======================================================
N_FARMS = 1600
N_TWO, N_THREE = 320, 40                     # farms surveyed twice / three times

iso_hh = rng.choice(list(COUNTRIES), size=N_FARMS, p=[COUNTRIES[k]["hh"] for k in COUNTRIES])
farm_ids = rng.choice(np.arange(100_000, 999_999), size=N_FARMS, replace=False)
farms = pd.DataFrame({
    "forestgardenid": [f"FG-{i:06d}" for i in farm_ids],
    "iso_trees": iso_hh,
})
farms["country"] = farms.iso_trees.map({k: v["name"] for k, v in COUNTRIES.items()})
# assign each farm to a project unit of its country
units_by_iso = {iso: UNITS[UNITS.iso == iso].reset_index(drop=True) for iso in COUNTRIES}
pick = [units_by_iso[iso].iloc[rng.integers(len(units_by_iso[iso]))] for iso in farms.iso_trees]
farms["project"] = [p.projectunit for p in pick]
farms["proj_area"] = [p.proj_area for p in pick]
farms["cell_id"] = [p.cell_id for p in pick]
farms["projectauto"] = farms["project"]
farms["village"] = [f"Village {rng.integers(1, 400):03d}" for _ in range(N_FARMS)]
# project (cohort) start date of the farm
farms["year_start_proj"] = rng.choice([2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023],
                                      size=N_FARMS, p=[.06, .03, .14, .07, .15, .25, .13, .10, .07])
farms["month_start_proj"] = rng.integers(1, 13, size=N_FARMS)
# latent household "wealth" and other farm-level characteristics
farms["u"] = rng.normal(0, 1, N_FARMS)
farms["n_surveys"] = 1
farms.loc[rng.choice(N_FARMS, N_TWO + N_THREE, replace=False), "n_surveys"] = 2
farms.loc[rng.choice(farms.index[farms.n_surveys == 2], N_THREE, replace=False), "n_surveys"] = 3

# survey timing (days since project start)
rows = []
for f in farms.itertuples():
    if f.n_surveys == 1:
        t = [rng.uniform(-30, 1400)]
    else:
        t = [rng.uniform(-30, 1000)]
        t.append(min(t[-1] + rng.uniform(300, 800), 1427))        # ~90% at least 365 days later
        if f.n_surveys == 3:
            t.append(min(t[-1] + rng.uniform(300, 500), 1427))
    for k, tt in enumerate(t):
        rows.append(dict(forestgardenid=f.forestgardenid, time_from_start=float(round(tt)), survey_k=k))
hh = pd.DataFrame(rows).merge(farms, on="forestgardenid")
n = len(hh)

# calendar dates of the surveys
start = pd.to_datetime(dict(year=hh.year_start_proj, month=hh.month_start_proj, day=1))
date = start + pd.to_timedelta(hh.time_from_start, unit="D")
hh["year_survey"] = date.dt.year.astype(float)
hh["month_survey"] = date.dt.month.astype(float)
hh["treat_year"] = np.clip(np.floor(np.maximum(hh.time_from_start, 0) / 365), 0, 3)
hh["post_dgp"] = (hh.survey_k > 0).astype(int)                  # within-farm "post first survey"

# --- household characteristics (farm-level, constant across surveys) -------
hh["participantage"] = np.clip(rng.normal(46, 12, n).round(), 18, 85)
hh.loc[rng.random(n) < 0.10, "participantage"] = np.nan
hh["participantgender"] = np.where(rng.random(n) < 0.57, "Male", "Female")
hh.loc[rng.random(n) < 0.05, "participantgender"] = ""
hh["gender"] = np.where(hh.participantgender == "", np.nan, (hh.participantgender == "Male").astype(float))
hh["typeoffarm"] = np.where(rng.random(n) < 0.98, "Individual", "Group Farm")
hh.loc[rng.random(n) < 0.10, "typeoffarm"] = ""
hh["individual"] = np.where(hh.typeoffarm == "", np.nan, (hh.typeoffarm == "Individual").astype(float))
hh["householdfemales"] = rng.poisson(4.3, n).astype(float)
hh["householdmales"] = rng.poisson(4.1, n).astype(float)
miss = rng.random(n) < 0.10
hh.loc[miss, ["householdfemales", "householdmales"]] = np.nan
hh["npeople"] = (hh.householdfemales + hh.householdmales).fillna(0)
hh.loc[hh.npeople == 0, "npeople"] = 1.0
hh["females_tot"] = np.where(hh.householdfemales.isna(), np.nan, hh.householdfemales / hh.npeople * 100)
hh["sizeofhousehold"] = np.where(rng.random(n) < 0.10, rng.integers(1, 20, n), np.nan)
hh["totalfarmland"] = np.where(rng.random(n) < 0.10, rng.uniform(0.2, 5, n).round(1), np.nan)
hh["income"] = np.where(rng.random(n) < 0.07, rng.normal(3e5, 1e6, n).round(), np.nan)
hh["farmingmainsource"] = yes_no(0.96, n, missing=0.10)
hh["guestsstaylastnight"] = yes_no(0.23, n, missing=0.10)
hh["householdchildren"] = yes_no(0.87, n, missing=0.10)
hh["childrenunder5"] = rng.choice(["none", "1-2 children", "3-5 children", "6-8 children", "9-12 children", "over 12 children", ""],
                                  size=n, p=[.17, .44, .15, .013, .004, .003, .22])
hh["childrennumber"] = rng.choice(["1-2 children", "3-5 children", "6-8 children", "9-12 children", "over 12 children", ""],
                                  size=n, p=[.24, .36, .13, .04, .01, .22])
CHILD_MAP = {"none": 0, "1-2 children": 1.5, "3-5 children": 4.5, "6-8 children": 7, "9-12 children": 10.5, "over 12 children": 13, "": np.nan}
hh["n_childrenunder5"] = hh.childrenunder5.map(CHILD_MAP)
hh["n_childrennumber"] = hh.childrennumber.map(CHILD_MAP)
hh["householdhighestschoollevel"] = rng.choice(["never attended", "some primary", "primary", "some secondary/tech", "secondary/tech", "higher", "do not know", ""],
                                               size=n, p=[.045, .105, .28, .18, .21, .08, .01, .09])
hh["childreninschool"] = yes_no(0.81, n, missing=0.22)
hh["mosteducated"] = np.where(rng.random(n) < 0.78, "Adult (older than 15)", "Child (15 or younger)")
hh.loc[rng.random(n) < 0.11, "mosteducated"] = ""
hh["adult_mostedu"] = np.where(hh.mosteducated == "", np.nan, (hh.mosteducated == "Adult (older than 15)").astype(float))
hh["reliantonotherhousehold"] = yes_no(0.37, n, missing=0.67)
hh["reliantothers"] = np.where(hh.reliantonotherhousehold == "", np.nan, (hh.reliantonotherhousehold == "Yes").astype(float))

# --- outcomes: livestock (TREATMENT EFFECT built in: post raises ownership by
#     ~6 pp and herd size by ~25 log points) ---------------------------------
u = hh.u.to_numpy()
post = hh.post_dgp.to_numpy()
p_own = 1 / (1 + np.exp(-(1.75 + 0.9 * u + 0.55 * post)))
hh["havelivestock"] = (rng.random(n) < p_own).astype(float)
hh["ownlivestock"] = np.where(hh.havelivestock == 1, "Yes", "No")
BINS = np.array([0, 2, 5.5, 9, 13, 18, 23.5, 26])
CUTS = np.array([1, 3.5, 7, 11, 15.5, 20.5, 25])
SPECIES = {"chicken": 6.5, "cattle": 1.9, "goat": 1.9, "sheep": 1.0, "mule": 0.55}
for s, lam in SPECIES.items():
    x = np.exp(np.log(lam) + 0.6 * u + 0.25 * post + rng.normal(0, 0.9, n))
    cnt = BINS[np.searchsorted(CUTS, x)]
    cnt = np.where(hh.havelivestock == 1, cnt, 0.0)
    hh[f"n_{s}count"] = cnt
    hh[f"{s}count"] = np.where(cnt == 0, "0", np.where(cnt == 2, "1-3", np.where(cnt == 5.5, "4-7", "8+")))
hh["tot_animal"] = hh[[f"n_{s}count" for s in SPECIES]].sum(axis=1)
PRICES = {"KEN": (5, 619, 116, 154, 46), "TZA": (10, 440, 220, 52, 48),
          "UGA": (10, 279, 935, 58, 52), "SEN": (8.3, 446, 424, 88, 49)}
pr = np.array([PRICES[i] for i in hh.iso_trees])
hh["animals_value"] = (pr[:, 0] * hh.n_chickencount + pr[:, 1] * hh.n_cattlecount + pr[:, 2] * hh.n_mulecount
                       + pr[:, 3] * hh.n_goatcount + pr[:, 4] * hh.n_sheepcount)

# --- farm size, trees, crops (TREATMENT EFFECT: trees +50 log points, crops +15) ---
size0 = np.exp(rng.normal(np.log(0.4), 0.6, N_FARMS))
size0[rng.random(N_FARMS) < 0.05] = 0.0                            # non-response coded as zero
size_map = dict(zip(farms.forestgardenid, size0))
hh["size"] = hh.forestgardenid.map(size_map).round(2)
hh["baselinesize"] = hh["size"]
hh.loc[rng.random(n) < 0.20, "baselinesize"] = np.nan
trees = np.exp(np.log(400) + 0.5 * u + 0.5 * post + rng.normal(0, 1.3, n)).round()
trees[rng.random(n) < 0.20] = 0
hh["trees"] = trees
hh.loc[rng.random(n) < 0.15, "trees"] = np.nan
hh["foodcrops"] = rng.poisson(5 * np.exp(0.2 * u + 0.15 * post)).astype(float)
hh["marketableproducts"] = rng.poisson(5.5 * np.exp(0.2 * u + 0.15 * post)).astype(float)
hh.loc[rng.random(n) < 0.15, ["foodcrops", "marketableproducts"]] = np.nan
hh["lnfoodcrops"] = np.log1p(hh.foodcrops)
hh["lnmarketableproducts"] = np.log1p(hh.marketableproducts)

# --- food security (TREATMENT EFFECT: small improvements) -----------------
hh["hdds"] = rng.binomial(12, np.clip(0.55 + 0.05 * u + 0.03 * post, 0.05, 0.95)).astype(float)
hh["hdds_cat"] = np.select([hh.hdds >= 8, hh.hdds >= 5], [1.0, 2.0], 3.0)
hh["hddscategory"] = hh.hdds_cat.map({1.0: "1. High Dietary Diversity", 2.0: "2. Moderate Dietary Diversity", 3.0: "3. Low Dietary Diversity"})
hh["hdds_id3"] = (hh.hdds_cat == 3).astype(np.int8)
hh["hfias"] = np.clip(rng.poisson(np.exp(np.log(6.8) - 0.25 * u - 0.15 * post)), 0, 27).astype(float)
hfcat = np.select([hh.hfias <= 1, hh.hfias <= 6, hh.hfias <= 13], [1, 2, 3], 4)
hh["hfiascategory"] = pd.Series(hfcat).map({1: "1. Food Secure", 2: "2. Mildly Food Insecure", 3: "3. Moderately Food Insecure", 4: "4. Severely Food Insecure"}).to_numpy()
hh["hfiascatid4"] = (hfcat == 4).astype(np.int8)
hh["hrcs_cat"] = rng.choice([1.0, 2.0, 3.0], size=n, p=[.04, .16, .80])
hh["hrcscategory"] = hh.hrcs_cat.map({1.0: "1. High Resilience Capacity", 2.0: "2. Moderate Resiliency Capacity", 3.0: "3. Low Resiliency Capacity"})
hh["worriedaboutfood"] = yes_no(0.63, n, missing=0.002)
hh["lackedfood"] = yes_no(0.67, n, missing=0.002)
hh["worriedfood"] = np.where(hh.worriedaboutfood == "", np.nan, (hh.worriedaboutfood == "Yes").astype(float))
hh["lackfood"] = np.where(hh.lackedfood == "", np.nan, (hh.lackedfood == "Yes").astype(float))
p_nofood = 1 / (1 + np.exp(-(-1.3 - 0.4 * u - 0.25 * post)))
hh["nofood"] = np.where(rng.random(n) < p_nofood, "Yes", "No")
hh["nofoodid"] = (hh.nofood == "Yes").astype(float)
hh["fewermeals"] = yes_no(0.49, n, missing=0.002)
hh["slepthungry"] = yes_no(0.21, n, missing=0.002)
hh["noeating24hrs"] = yes_no(0.14, n, missing=0.002)

# --- local GDP (cell x survey year) ---------------------------------------
hh = hh.merge(CELLS.rename(columns={"year": "year_survey"}), on=["cell_id", "year_survey"], how="left")
for v in ["gcpc_2017usd", "ch_gcpc_2017usd"]:
    hh[f"w1_{v}" if v == "gcpc_2017usd" else "w1ch_gcpc_2017usd"] = winsor(hh[v], 1, 99)
    hh[f"w5_{v}" if v == "gcpc_2017usd" else "w5ch_gcpc_2017usd"] = winsor(hh[v], 5, 95)
hh["pop_cell"] = hh.pop_cell

# --- winsorised versions and zero-fills, as in 3_put_files_together ------
hh["w_trees"] = winsor(hh.trees.fillna(0), 1, 99)
hh["w_marketableproducts"] = winsor(hh.marketableproducts.fillna(0), 1, 99)
hh["w_foodcrops"] = winsor(hh.foodcrops.fillna(0), 1, 99)
hh["w_animals_value"] = winsor(hh.animals_value, 1, 99)
hh["w_lntrees"] = winsor(np.log1p(hh.trees), 1, 99)
hh["w_npeople"] = winsor(hh.npeople, 1, 99)
hh["w_size"] = winsor(hh["size"], 1, 99)
hh["w_participantage"] = winsor(hh.participantage, 1, 99)
hh["w_lnbaselinetrees"] = winsor(np.log1p(hh.trees), 1, 99)

# per-hectare variables (missing when size==0, as division by zero in Stata)
sz = hh["size"].replace(0, np.nan)
hh["trees_per_ha"] = hh.trees / sz
hh["mark_prod_ha"] = hh.marketableproducts / sz
hh["foodcrop_ha"] = hh.w_foodcrops / sz
for v in ["trees_per_ha", "mark_prod_ha", "foodcrop_ha"]:
    hh[f"w1_{v}"] = winsor(hh[v], 1, 99)
    hh[f"lnw1_{v}"] = np.log1p(hh[f"w1_{v}"])
hh["tot_animal_ha"] = hh.tot_animal / sz
hh["wtot_animal_ha"] = winsor(hh.tot_animal_ha, 1, 99)
hh["lnwtot_animal_ha"] = np.log1p(hh.wtot_animal_ha)
hh["tot_animal_w"] = winsor(hh.tot_animal, 1, 99)

# --- baseline (first-survey) controls, wb_* ---------------------------------
hh = hh.sort_values(["forestgardenid", "time_from_start"]).reset_index(drop=True)
first = hh.groupby("forestgardenid").head(1).set_index("forestgardenid")
for src, dst in [("gender", "wb_gender"), ("adult_mostedu", "wb_adult_mostedu"), ("w_npeople", "wb_w_npeople"),
                 ("females_tot", "wb_females_tot"), ("w_size", "wb_w_size"), ("worriedfood", "wb_worriedfood"),
                 ("lackfood", "wb_lackfood"), ("w_participantage", "wb_w_participantage"), ("pop_cell", "wb_pop_cell"),
                 ("w1_gcpc_2017usd", "b_w1_gcpc_2017usd"), ("n_childrennumber", "b_n_childrennumber"),
                 ("n_childrenunder5", "b_n_childrenunder5"), ("hdds", "wb_hdds"), ("trees", "wb_trees")]:
    hh[dst] = hh.forestgardenid.map(first[src])
# script 3 zero-fills these controls where missing
for v in ["havelivestock", "tot_animal", "w_trees", "w_marketableproducts", "w_foodcrops", "nofoodid", "hfiascatid4",
          "hdds", "size", "npeople", "wb_adult_mostedu", "wb_gender", "wb_w_npeople", "wb_w_size", "wb_w_participantage"]:
    hh[v] = hh[v].fillna(0)

# --- fixed-effect group ids (egen group()) ----------------------------------
hh["fe_farm"] = group_id(hh, ["forestgardenid"]).astype(float)
hh["fe_proj"] = group_id(hh, ["project"]).astype(float)
hh["fe_proj_area"] = group_id(hh, ["proj_area"]).astype(float)
hh["fe_cell"] = group_id(hh, ["cell_id"]).astype(float)
hh["fe_projyr"] = group_id(hh, ["project", "year_survey"]).astype(float)
hh["fe_country"] = group_id(hh, ["country"]).astype(float)
hh["fe_ctryyr"] = group_id(hh, ["country", "year_survey"]).astype(float)
hh["fe_ctryyrmo"] = group_id(hh, ["country", "year_survey", "month_survey"]).astype(float)
hh["fe_moyr"] = group_id(hh, ["year_survey", "month_survey"]).astype(float)
hh["fe_ctrymo"] = group_id(hh, ["country", "month_survey"]).astype(float)
hh["fe_year"] = group_id(hh, ["year_survey"]).astype(float)
hh["post1"] = (hh.treat_year > 1).astype(float)
hh["year_cohort"] = hh.year_survey - hh.treat_year

# --- sample-definition flags, EXACTLY as in 3_put_files_together_v11.do ----
g = hh.groupby("fe_farm")
hh["n_obs"] = g.fe_farm.transform("size").astype(float)
hh["first_obs"] = (g.cumcount() == 0).astype(np.int8)
hh["ref_time"] = g.time_from_start.transform("first")
hh["gap_from_ref"] = hh.time_from_start - hh.ref_time
hh["gap_ok"] = ((hh.first_obs == 1) | ((hh.first_obs == 0) & (hh.gap_from_ref >= 365))).astype(np.int8)
hh["last_obs"] = (g.cumcount(ascending=False) == 0).astype(np.int8)
hh["gap_ok_fl"] = (hh.gap_ok.astype(bool) & ((hh.first_obs == 1) | (hh.last_obs == 1))).astype(np.int8)
hh["post"] = (hh.first_obs == 0).astype(np.int8)
hh["n_post"] = g.post.transform("sum").astype(float)
hh["n_pre"] = g.first_obs.transform("sum").astype(float)
hh["insample"] = ((hh.n_post >= 1) & (hh.n_pre >= 1)).astype(np.int8)

# first-to-last changes (missing for single-survey farms)
def first_last_change(v, relative=False):
    f_ = g[v].transform("first"); l_ = g[v].transform("last")
    out = np.where(hh.n_obs > 1, (l_ - f_) / f_ if relative else l_ - f_, np.nan)
    return out

for v in ["tot_animal", "tot_animal_ha", "foodcrop_ha", "mark_prod_ha", "trees_per_ha", "hfiascatid4", "nofoodid", "hdds"]:
    hh[f"ch_{v}"] = first_last_change(v)
    hh[f"wch_{v}"] = winsor(hh[f"ch_{v}"], 1, 99)
for v in ["w1_gcpc_2017usd", "w5_gcpc_2017usd"]:
    hh[f"ch_{v}"] = first_last_change(v, relative=True)
hh["wch_trees1000_per_ha"] = hh.wch_trees_per_ha / 1000
# GCP per capita expressed in USD (the processed file multiplies by 1e9)
hh["w1_gcpc_2017usd"] = hh.w1_gcpc_2017usd * 1e9
hh["w5_gcpc_2017usd"] = hh.w5_gcpc_2017usd * 1e9

LABELS_HH = {
    "forestgardenid": "Forest Garden", "projectauto": "Project (auto)",
    "householdfemales": "Household Females", "hdds": "HDDS", "hfias": "HFIAS",
    "hfiascategory": "HFIAS Category", "tot_animal": "N. Livestock", "havelivestock": "1(Own livestock)",
    "size": "Farm size (ha)", "w_trees": "N. trees",
    "w1_gcpc_2017usd": "GCP per Capita", "w5_gcpc_2017usd": "GCP per Capita",
    "w1ch_gcpc_2017usd": "Change in GCP per Capita", "w5ch_gcpc_2017usd": "Change in GCP per Capita",
    "wb_gender": "1(House. Head Male)", "wb_adult_mostedu": "1(Adult Most Edu.)", "wb_w_npeople": "N. House. Members",
    "wb_w_size": "Farm Size (ha)", "wb_w_participantage": "House. Head Age", "w1_trees_per_ha": "N. Trees per ha",
    "w1_mark_prod_ha": "N. Marketable Crops per ha", "w1_foodcrop_ha": "N. Food Crops per ha",
    "wtot_animal_ha": "tot_animal_ha", "wch_trees1000_per_ha": "Ch. Trees per ha (1000s/ha)",
    "wch_mark_prod_ha": "Ch. Market Crops per ha", "wch_foodcrop_ha": "Ch. Food Crops per ha", "post": "1(Post 1st Survey)",
}

# --- variables written to the file: ONLY those used by the paper's scripts ---
# (02_summary_stats.do, 03_main_estimates.do, 06_python_replication.py). Every
# other variable built above is an intermediate of the data-generating process
# (or a raw-survey field the paper does not use) and is not written. Variables
# that a script only labels, e.g. `gender`, are kept so the scripts run unchanged.
KEEP_HH = [
    # identifiers and timing
    "forestgardenid", "iso_trees", "project", "proj_area", "projectauto", "time_from_start",
    "year_survey", "month_survey", "treat_year", "year_cohort",
    # fixed-effect ids (egen group())
    "fe_farm", "fe_proj_area", "fe_projyr", "fe_country", "fe_year",
    # sample-definition flags (3_put_files_together)
    "n_obs", "first_obs", "last_obs", "gap_ok", "gap_ok_fl", "post", "insample",
    # household characteristics (summary statistics, World Bank and single/multi comparisons)
    "gender", "householdfemales", "npeople", "females_tot", "farmingmainsource", "guestsstaylastnight",
    "householdchildren", "n_childrenunder5", "n_childrennumber", "householdhighestschoollevel",
    "childreninschool", "reliantothers",
    # livestock outcomes
    "havelivestock", "n_chickencount", "n_cattlecount", "n_goatcount", "n_sheepcount", "n_mulecount",
    "tot_animal", "tot_animal_ha", "wtot_animal_ha", "w_animals_value",
    # farm size, trees and crops
    "size", "baselinesize", "trees", "foodcrops", "marketableproducts", "w_trees", "w_marketableproducts",
    "w_foodcrops", "trees_per_ha", "w1_trees_per_ha", "w1_mark_prod_ha", "w1_foodcrop_ha",
    # food security
    "hdds", "hdds_id3", "hfias", "hfiascategory", "hfiascatid4", "worriedfood", "nofoodid", "fewermeals",
    "slepthungry", "noeating24hrs",
    # local GDP (gross cell product) controls
    "gcp_2017usd", "l1_gcpc_2017usd", "w1_gcpc_2017usd", "w5_gcpc_2017usd", "w1ch_gcpc_2017usd",
    "w5ch_gcpc_2017usd", "ch_w5_gcpc_2017usd",
    # baseline (first-survey) controls
    "wb_gender", "wb_adult_mostedu", "wb_w_npeople", "wb_females_tot", "wb_w_size", "wb_worriedfood",
    "wb_lackfood", "wb_w_participantage", "b_n_childrennumber",
    # first-to-last changes (median regressions, channel table)
    "ch_tot_animal", "wch_tot_animal_ha", "wch_foodcrop_ha", "wch_mark_prod_ha", "wch_trees_per_ha",
    "wch_trees1000_per_ha",
]
assert len(KEEP_HH) == len(set(KEEP_HH)) == 88, "KEEP_HH must list 88 distinct variables"
missing = [v for v in KEEP_HH if v not in hh.columns]
assert not missing, f"KEEP_HH names not built above: {missing}"
hh_out = hh[KEEP_HH]
hh_out.to_stata(FINAL / "TREES_data_pseudo.dta", write_index=False, version=118,
                variable_labels={k: v for k, v in LABELS_HH.items() if k in hh_out})
print(f"TREES_data_pseudo: {len(hh_out):,} rows, {hh_out.forestgardenid.nunique():,} farms, "
      f"{hh_out.shape[1]} vars (of {hh.shape[1]} built); estimation sample rows = "
      f"{int(((hh_out.insample==1)&(hh_out.n_obs>1)&(hh_out.gap_ok==1)&(hh_out.gap_ok_fl==1)).sum())}")


# =======================================================
# 2. Stand-in for the World Bank LSMS-ISA Tanzania benchmark (final/tanzania_merged.dta)
#    The real file is built from the public LSMS-ISA wave 5 microdata; this is a
#    random stand-in with the same columns so wb_comparison runs.
# =======================================================
m = 1000
tz = pd.DataFrame({
    "region": rng.choice(["manyara", "singida", "tabora", "dodoma", "mwanza", "arusha"], size=m, p=[.2, .18, .17, .15, .15, .15]),
    "npeople": rng.poisson(5, m) + 1.0, "year_survey": rng.choice([2020.0, 2021.0], size=m), "month_survey": rng.integers(1, 13, m).astype(float),
})
tz["householdfemales"] = np.clip(rng.binomial(tz.npeople.astype(int), 0.5), 0, None).astype(float)
tz["p_females"] = tz.householdfemales / tz.npeople
tz["n_childrennumber"] = rng.poisson(2.5, m).astype(float)
tz["n_childrenunder5"] = rng.poisson(0.8, m).astype(float)
tz["p_child5yrs"] = tz.n_childrenunder5 / tz.npeople
tz["farm_area"] = np.exp(rng.normal(np.log(1.5), 0.8, m)).round(2)
tz["trees"] = np.where(rng.random(m) < 0.6, np.exp(rng.normal(np.log(30), 1.5, m)).round(), 0.0)
tz["trees_per_ha"] = tz.trees / tz.farm_area
tz["havelivestock"] = (rng.random(m) < 0.62).astype(float)
tz["tot_animal"] = np.where(tz.havelivestock == 1, rng.poisson(9, m), 0).astype(float)
tz["tot_animal_ha"] = tz.tot_animal / tz.farm_area
tz["n_crops"] = rng.poisson(3, m).astype(float)
tz["hfias"] = np.clip(rng.poisson(5, m), 0, 27).astype(float)
tz["hfias_category"] = np.select([tz.hfias <= 1, tz.hfias <= 6, tz.hfias <= 13], [0.0, 1.0, 2.0], 3.0)
tz["hfiascatid4"] = (tz.hfias_category == 3).astype(float)
for v, p in [("worriedfood", .55), ("reduce_n_meals", .35), ("nofoodid", .18), ("farm_main_income", .8), ("child_in_school", .85)]:
    tz[v] = (rng.random(m) < p).astype(float)
tz["highest_sch_level"] = rng.choice([0.0, 1.0, 2.0, 3.0], size=m, p=[.15, .55, .25, .05])
tz.to_stata(FINAL / "tanzania_merged.dta", write_index=False, version=118)
print(f"tanzania_merged: {len(tz):,} rows")

# =======================================================
# 3. Satellite (EVI) farm-year panel
# =======================================================
EVI_YEARS = list(range(2017, 2025))
COHORT_P = {2018: .05, 2019: .08, 2020: .18, 2021: .10, 2022: .07, 2023: .27, 2024: .25}
# dynamic treatment effect on EVI (raw 0-1 scale) by event time (no pre-trend;
# the scripts keep event times -3..+4)
TAU = {-3: 0.0, -2: 0.0, -1: 0.0, 0: 0.003, 1: 0.010, 2: 0.015, 3: 0.020, 4: 0.025, 5: 0.025}


def sample_points(iso, k, lat0, lon0):
    """fictional coordinates: scatter around the (fictional) project-area centre"""
    return lat0 + rng.normal(0, 0.12, k), lon0 + rng.normal(0, 0.12, k)


def make_evi_panel(n_treated, n_control, seed):
    r = np.random.default_rng(seed)
    iso_p = [COUNTRIES[k]["sat"] for k in COUNTRIES]
    # --- TREES farms (FG-): all treated (the paper no longer uses the
    #     signed-up-but-never-treated farms, so none are generated) ---
    n_fg = n_treated
    fg_ids = r.choice(np.arange(100, 999_999), size=n_fg, replace=False)
    fg = pd.DataFrame({"forestgardenid": [f"FG-{i:06d}" for i in fg_ids],
                       "iso_trees": r.choice(list(COUNTRIES), size=n_fg, p=iso_p)})
    pick = [units_by_iso[iso].iloc[r.integers(len(units_by_iso[iso]))] for iso in fg.iso_trees]
    fg["projectunit"] = [p.projectunit for p in pick]
    fg["proj_area"] = [p.proj_area for p in pick]
    fg["cell_id"] = [p.cell_id for p in pick]
    fg["lat"], fg["lon"] = zip(*[(a, b) for a, b in zip(*sample_points(None, n_fg, np.array([p.lat0 for p in pick]), np.array([p.lon0 for p in pick])))])
    fg["TREES_sample"] = 1.0
    fg["type"] = ""
    fg["original_farm_id"] = ""
    fg["status"] = r.choice(["Open", "Closed: Graduated"], size=n_fg, p=[.77, .23])
    fg["calyear_treat"] = r.choice(list(COHORT_P), size=n_fg, p=list(COHORT_P.values())).astype(float)
    # --- random control farms (CF-), each seeded from a TREES farm ---
    seeds = fg.sample(n_control, replace=True, random_state=int(r.integers(1e9))).reset_index(drop=True)
    cf = seeds[["iso_trees", "projectunit", "proj_area", "cell_id"]].copy()
    cf["forestgardenid"] = [f"CF-{i}" for i in r.choice(np.arange(1000, 99_999), size=n_control, replace=False)]
    cf["original_farm_id"] = seeds.forestgardenid
    cf["lat"] = np.nan; cf["lon"] = np.nan
    cf["TREES_sample"] = 0.0; cf["type"] = "control"; cf["status"] = "Random Control Farms"; cf["calyear_treat"] = 0.0
    allf = pd.concat([fg, cf], ignore_index=True)
    allf["trees_treated_farms"] = (allf.calyear_treat > 0).astype(float)
    allf["never_treat"] = 1.0 - allf.trees_treated_farms
    allf["iso"] = allf.iso_trees
    # project-unit cohort: first TREES cohort in the unit (controls inherit it)
    unit_coh = allf[allf.calyear_treat > 0].groupby("projectunit").calyear_treat.min()
    allf["projunit_calyear_treat"] = allf.projectunit.map(unit_coh)
    allf["projunit_calyear_treat"] = allf.projunit_calyear_treat.fillna(allf.projectunit.map(UNITS.set_index("projectunit").unit_start).astype(float))
    # farm-level EVI heterogeneity
    allf["alpha"] = np.clip(r.normal(0.26, 0.09, len(allf)), 0.08, 0.60)
    # --- farm x year panel ---
    panel = allf.loc[allf.index.repeat(len(EVI_YEARS))].reset_index(drop=True)
    panel["year"] = np.tile(EVI_YEARS, len(allf)).astype(np.int16)
    N = len(panel)
    cy = {(c, y): r.normal(0, 0.02) for c in COUNTRIES for y in EVI_YEARS}
    panel["treat_year"] = np.where(panel.calyear_treat > 0, panel.year - panel.calyear_treat, np.nan)
    tau = panel.treat_year.map(TAU).fillna(0.0).to_numpy()
    panel = panel.merge(CELLS, on=["cell_id", "year"], how="left")
    evi = (panel.alpha.to_numpy() + np.array([cy[(c, y)] for c, y in zip(panel.iso_trees, panel.year)])
           + tau + 0.002 * (panel.gcpc_2017usd.to_numpy() * 1e6 - 1.5) + r.normal(0, 0.03, N))
    panel["evi_median"] = np.clip(evi, -0.06, 0.75)
    panel["evi_mean"] = panel.evi_median + r.normal(0.02, 0.03, N)
    blow = r.random(N) < 0.0005                                   # a few blown-up pixels, as in the raw extraction
    panel.loc[blow, "evi_mean"] = r.choice([-1, 1], size=blow.sum()) * 10 ** r.uniform(3, 9, blow.sum())
    # balanced panel: every farm has all 8 years (the scripts apply the -3..+4
    # event window themselves), so the file has exactly n_farms x 8 rows
    panel["n_obs"] = 8.0
    # GDP controls (winsorised), EVI in x100 units
    panel["gcpc_2017usd_w1"] = winsor(panel.gcpc_2017usd, 1, 99)
    panel["gcpc_2017usd_w5"] = winsor(panel.gcpc_2017usd, 5, 95)
    panel["ch_gcpc_2017usd_w1"] = winsor(panel.ch_gcpc_2017usd, 1, 99)
    panel["w1_evi_median"] = winsor(panel.evi_median, 1, 99) * 100
    panel["w5_evi_median"] = winsor(panel.evi_median, 5, 95)
    panel["w1_evi_mean"] = winsor(panel.evi_mean, 1, 99) * 100
    p99, p1 = np.nanpercentile(panel.evi_median, [99, 1])
    panel["to_drop"] = panel.groupby("forestgardenid").evi_median.transform(lambda s: float(((s >= p99) | (s <= p1)).any()))
    for k in range(5, 14):
        panel[f"treat_yearid{k}"] = (panel.treat_year == k - 8).astype(np.int8)
    panel["treated"] = ((panel.calyear_treat > 0) & (panel.year > panel.calyear_treat)).astype(float)
    panel["post"] = panel.treated
    panel["treat_post"] = panel.treated
    panel["calyear"] = panel.year.astype(float)
    panel["fe_proj"] = group_id(panel, ["projectunit"]).astype(float)
    panel["fe_proj_area"] = group_id(panel, ["proj_area"]).astype(float)
    panel["fe_cell"] = group_id(panel, ["cell_id"]).astype(float)
    panel["fe_country"] = group_id(panel, ["iso_trees"]).astype(float)
    panel["fe_ctryyr"] = group_id(panel, ["iso_trees", "year"]).astype(float)
    panel["fe_farm"] = group_id(panel, ["forestgardenid"]).astype(float)
    panel["fe_year"] = group_id(panel, ["year"]).astype(float)
    panel = panel.drop(columns=["alpha"]).sort_values(["forestgardenid", "year"]).reset_index(drop=True)
    return panel


EVI_LABELS = {"w1_evi_median": "EVI median (x 100)", "w1_evi_mean": "EVI mean (x 100)", "gcpc_2017usd": "GCP per Capita",
              "ch_gcpc_2017usd": "Change in GCP per Capita", "calyear_treat": "Calendar year of treatment (0 = never)"}
ev = make_evi_panel(n_treated=2600, n_control=2400, seed=SEED + 11)   # 5,000 farms x 8 years
assert len(ev) == 40_000 == ev.forestgardenid.nunique() * len(EVI_YEARS), "EVI panel must be balanced, 40,000 rows"
ev.to_stata(FINAL / "EVI_farms_pseudo.dta", write_index=False, version=118,
            variable_labels={k: v for k, v in EVI_LABELS.items() if k in ev})
print(f"EVI_farms_pseudo: {len(ev):,} farm-years, {ev.forestgardenid.nunique():,} farms "
      f"({int(ev.drop_duplicates('forestgardenid').trees_treated_farms.sum()):,} treated), {ev.shape[1]} vars")

print("\nDone. Pseudo data written to", FINAL)
