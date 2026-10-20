# Code and demo data for *Natural capital restoration and household wealth in sub-Saharan Africa*

Claudio Rizzi (IESE Business School), Simon Xu (Stanford University), Geoffrey Heal (Columbia University).
Corresponding authors: Claudio Rizzi (crizzi@iese.edu) and Simon Xu (sdzxu@stanford.edu).

This package contains the analysis code of the manuscript and a **synthetic demo
dataset** on which every script runs end-to-end. It follows the structure asked
for in Nature's *Code and Software Submission Checklist*: system requirements,
installation guide, demo (with expected output and run time), instructions for
use, and reproduction instructions.

The package (code, synthetic demo data and the expected output of the demo) is
publicly available at https://github.com/claudiorizzi/natural_capital_restoration.

---------------------------------------------------------------------------------

## Contents

```
natural_capital_restoration/   repository root (the folder holding this README)
├── README.md                  this file
├── LICENSE                    MIT (code); CC0 (synthetic data)
├── pyproject.toml, uv.lock    Python environment (pinned versions)
├── .env                       paths used by the Python scripts (relative, ready to use)
├── code/
│   ├── run_demo.do                master do-file: runs 02 and 03 on the demo data
│   ├── 01_make_pseudo_data.py     generates the synthetic demo data (documents the data structure)
│   ├── 02_summary_stats.do        summary statistics and sample-comparison tables
│   ├── 03_main_estimates.do       all estimations: household outcomes (OLS / PPML with farm
│   │                              fixed effects), satellite EVI event studies (stacked DiD and
│   │                              Sun & Abraham), robustness checks
│   ├── 04_maps.py                 map of farm locations (hexbin over Natural Earth boundaries)
│   ├── 05_value_scc.py            present value of the climate benefits (social cost of carbon)
│   ├── 06_python_replication.py   numpy/pandas re-implementation of the core estimators (cross-check)
│   ├── gee_control_farms_50km.js  Google Earth Engine: random control locations (1-50 km); upstream
│   │                              of the data, not numbered (runs in the GEE code editor)
│   └── gee_evi_extraction.js      Google Earth Engine: annual EVI per farm (Sentinel-2, Cloud Score+)
├── data/
│   ├── final/                 synthetic demo data (3 Stata .dta files, about 17 MB)
│   └── country_boundaries/    Natural Earth 1:110m admin-0 boundaries (public domain)
└── output/                    written by the demo: tables/ (LaTeX), figures/ (PNG), logs
```

The numbered scripts run in the order of their prefix (01 to 06; the master
do-file runs 02 and 03). The two Google Earth Engine scripts are upstream of
the data and are not part of the runnable pipeline, so they carry no number.
The two Stata scripts contain exactly the analyses reported in the manuscript
and its Supplementary Information (every program is mapped to a figure or
table in the file headers and in section 3 below). The Python and Google Earth
Engine scripts are the ones used for the manuscript; their only edits are the
data and output paths, which are read from `.env`.

---------------------------------------------------------------------------------

## Data availability and the synthetic demo data

The household-survey data and the farm-location data were provided by Trees for
the Future (TREES) under a data-use agreement and **cannot be redistributed**:

* the agreement does not permit onward sharing;
* the records are personal data on individual rural households (household
  composition, food security, livestock, income sources);
* every farm carries GPS coordinates, and the satellite panel is built from
  them, so the data cannot be anonymised without destroying the information the
  analysis relies on (location is the key to the satellite outcomes and to the
  control-location draw).

Researchers may request the data directly from TREES (trees.org), subject to its
confidentiality requirements.

The manuscript also uses two public datasets:

* the **World Bank LSMS-ISA Tanzania National Panel Survey 2020-21 (Wave 5)**
  (reference ID TZA_2020_NPS-R5_v02_M), the household benchmark of
  Supp. Table 11. It can be downloaded, after free registration, from the World
  Bank Microdata Library: https://microdata.worldbank.org/index.php/catalog/5639;
* the **local (subnational) GDP grid of Rossi-Hansberg and Zhang**, *Local
  Economies, Global Perspective: Illuminating Subnational GDP Worldwide* (BFI
  Data Studio, University of Chicago), the source of the gross-cell-product per
  capita controls. The grids (ZIP files with shapefiles and CSV files at 1, 0.5
  and 0.25 degrees; the paper uses the 0.25-degree grid) can be downloaded from
  https://bfidatastudio.org/project/local-economies-global-perspective-illuminating-subnational-gdp-worldwide/.

For simplicity of the replication code, and because the farm data cannot be
shared, the package does not include any raw input files (neither the
restricted TREES files nor the two public datasets). The final datasets in
`data/final/` contain instead synthetic versions of the variables used in the
paper, including the variables derived from the two public sources.

To make the code runnable, `code/01_make_pseudo_data.py` generates **pseudo
datasets with the variable names, types, panel structure and sample-definition
flags of the processed files the scripts read**, limited to the variables the
scripts use. Every value is simulated; no farm, household, identifier, project
name or coordinate is real. The pseudo data are a synthetic version of the
data, not a copy: **the numbers of farms, surveys and farm-years do not match
the samples of the manuscript** (they are set for a quick demo). The
data-generating process builds in treatment effects so that the demo produces
sensible output, but **estimates obtained on the demo data are not the
manuscript's estimates**.

| demo file (data/final/)   | stands in for                                                    | size |
|---------------------------|------------------------------------------------------------------|------|
| `TREES_data_pseudo.dta`   | household survey panel (all household estimates, summary statistics, sample comparisons) | 2,000 farm-surveys, 1,600 farms (360 surveyed twice or more), 88 variables (only those the scripts use) |
| `EVI_farms_pseudo.dta`    | satellite farm-year panel (EVI event studies)                    | 40,000 farm-years: 5,000 farms (2,600 treated farms, 2,400 random control locations) x 8 years, 2017-2024 |
| `tanzania_merged.dta`     | World Bank LSMS-ISA Tanzania benchmark (public data; random stand-in) | 1,000 households |

---------------------------------------------------------------------------------

## 1. System requirements

**Operating system.** Tested on Windows 11 Enterprise (x64). Stata and Python are
cross-platform; all paths in the package are relative, so macOS and Linux should
work unchanged (not tested).

**Stata** (all estimations, tables and figures): tested with StataNow 19.5 MP;
any Stata 16 or later should work (`jwdid` requires 16+). User-written packages
from SSC, with the versions used:

| package            | version tested      | used for                                           |
|--------------------|---------------------|----------------------------------------------------|
| ftools             | 2.50.0 (09jan2026)  | dependency of reghdfe / ppmlhdfe                   |
| reghdfe            | 6.13.1 (10jan2026)  | OLS with high-dimensional fixed effects            |
| ppmlhdfe           | 2.3.3 (02nov2025)   | Poisson pseudo-maximum likelihood with fixed effects|
| jwdid              | 2.21                | Wooldridge (2021) extended TWFE ("stacked DiD")    |
| hdfe               | 3.2.9               | dependency of jwdid                                |
| eventstudyinteract | 0.1 (24jan2022)     | Sun & Abraham (2021) interaction-weighted estimator|
| avar               | 1.0.07              | dependency of eventstudyinteract                   |
| estout / esttab    | 2.1.1 (10jun2022)   | LaTeX tables                                       |
| erepost            | 1.0.3 (26jun2025)   | re-posting event-study coefficients                |
| distinct           | 1.2.1               | counts of farms                                    |
| winsor2            | 1.1 (2014)          | winsorisation                                      |
| qreg2              | (SSC, 2016)         | median regressions with clustered SE               |
| psmatch2           | 4.0.12 (30jan2016)  | matched-control robustness                         |

**Python** (maps, carbon valuation, demo data, cross-check): Python 3.11 or later
(tested 3.12.3) managed with [uv](https://docs.astral.sh/uv/) (tested uv 0.9.18).
Packages (pinned in `uv.lock`; tested versions): pandas 3.0.6, numpy 2.4.6,
scipy 1.17.1, geopandas 1.2.0, shapely 2.2.0, pyproj 3.7.2, matplotlib 3.11.2,
python-dotenv 1.2.4.

**Google Earth Engine** (only to re-extract the satellite data, not needed for
the demo): a GEE account; the two `.js` files are pasted into the GEE code
editor (code.earthengine.google.com).

**Hardware.** No non-standard hardware. The demo runs on any desktop or laptop
with 4 GB of RAM. The manuscript's satellite panel (about 300,000 farm-years)
needs roughly 2 GB of RAM in Stata; with the manuscript's 1,000 placebo
re-estimations the placebo-timing test is the slowest step (the package
ships 5 draws, see section 5).

---------------------------------------------------------------------------------

## 2. Installation guide

1. **Get the package** (about 17 MB):
   ```bash
   git clone https://github.com/claudiorizzi/natural_capital_restoration.git
   ```
   or download the ZIP from the repository page and unpack it. The folder
   `natural_capital_restoration` is the repository root referred to below.
2. **Stata packages** (once, from the Stata command line; about 2 minutes):
   ```stata
   ssc install ftools
   ssc install reghdfe
   ssc install hdfe
   ssc install ppmlhdfe
   ssc install jwdid
   ssc install avar
   ssc install eventstudyinteract
   ssc install estout
   ssc install erepost
   ssc install distinct
   ssc install winsor2
   ssc install qreg2
   ssc install psmatch2
   ```
3. **Python environment** (once; 1-2 minutes, about 350 MB download):
   install uv (https://docs.astral.sh/uv/getting-started/installation/), then
   from the repository root run
   ```bash
   uv sync
   ```
   This creates `.venv` with the pinned package versions. `uv run python ...`
   (used below) does this automatically on first use, so step 3 can be skipped.

Typical total install time on a normal desktop: under 5 minutes.

---------------------------------------------------------------------------------

## 3. Demo

All commands are run **from the repository root** (the folder holding this
README).

**Step 0 (optional): regenerate the demo data** (about 20 seconds; the shipped
files were produced with this command and the fixed seed 20261019):
```bash
uv run python code/01_make_pseudo_data.py
```

**Step 1: all Stata estimations.** Open Stata, `cd` to the repository root
and run
```stata
do code/run_demo.do
```
or, in batch mode, `"C:\Program Files\Stata19\StataMP-64.exe" /e do code\run_demo.do`
(Windows) / `stata-mp -b do code/run_demo.do` (Mac, Linux). The master file
defines the programs of `02_summary_stats.do` and `03_main_estimates.do` and
calls them one by one, in that order. Each program is wrapped in `capture noisily`, so a
failure in one does not stop the others; the return code and run time of every
program are printed at the end and saved to `output/run_demo_summary.csv`. The
full log is `output/run_demo.log`.

**Step 2: map of farm locations** (about 10 seconds):
```bash
uv run python code/04_maps.py
```

**Step 3: carbon valuation** (about 1 second; uses no survey data, so its
output is the manuscript's own table):
```bash
uv run python code/05_value_scc.py
```

**Step 4: independent cross-check of the Stata results** (about 5 seconds):
```bash
uv run python code/06_python_replication.py
```

### Expected run time (Step 1, demo data, StataNow 19.5 MP on a laptop)

| program | script | seconds |
|---|---|---:|
| sum_stats_farms | 02_summary_stats.do | 0.1 |
| sum_stats_evi | 02_summary_stats.do | 0.1 |
| wb_comparison | 02_summary_stats.do | 2.1 |
| single_multi_comparison | 02_summary_stats.do | 0.2 |
| baseline | 03_main_estimates.do | 3.6 |
| leave_one_out | 03_main_estimates.do | 6.0 |
| median_reg | 03_main_estimates.do | 0.4 |
| plots_wealth | 03_main_estimates.do | 1.4 |
| evi_analysis | 03_main_estimates.do | 32.2 |
| evi_robustness | 03_main_estimates.do | 34.9 |
| trees_crops | 03_main_estimates.do | 1.9 |
| channel | 03_main_estimates.do | 0.5 |
| health | 03_main_estimates.do | 1.2 |
| add_rob | 03_main_estimates.do | 1.3 |
| **total** | | **85.9** |

Times from `output/run_demo_summary.csv` of the shipped run (StataNow 19.5 MP,
Windows 11 laptop), with the package's 5 placebo draws; every program returned
code 0.

### Expected output

Step 1 writes LaTeX tables to `output/tables/` and PNG figures to
`output/figures/`. Every output corresponds to a figure or table of the
manuscript (Fig.) or of the Supplementary Information (Supp.):

| program (02_summary_stats.do) | output                                                                 | manuscript |
|-------------------------------|------------------------------------------------------------------------|------------|
| sum_stats_farms, sum_stats_evi | summary_stats.tex, summary_stats_EVI.tex                    | Supp. Table 3: summary statistics |
| wb_comparison       | comparison_table.tex                                                   | Supp. Table 11: TREES farms vs World Bank LSMS-ISA households |
| single_multi_comparison | comparison_single_multi.tex                                        | Supp. Table 12: single- vs multi-survey farms |
| **(03_main_estimates.do)** |                                                                    |            |
| baseline            | own_n_livestock.tex                                                    | Supp. Table 4: livestock ownership and herd size, within-farm OLS / PPML |
|                     | rob_alternative_fes.tex                                                | Supp. Table 8: alternative fixed effects, single-survey farms, livestock per hectare, herd value |
| leave_one_out       | leave_one_out_animal.png, leave_one_out.png                            | Supp. Figs 5 and 4: leave-one-animal-out, leave-one-country-out |
| median_reg          | rob_median_reg.tex                                                     | Supp. Table 10: median regressions |
| plots_wealth        | wealth_post_coefficients.png                                           | Fig. 2: post coefficients across wealth outcomes |
| evi_analysis        | evi_median_mean.tex, evi_median_combined.png, evi_mean_combined.png    | Supp. Table 7 and Fig. 1: EVI event studies, stacked DiD (jwdid) and Sun & Abraham |
| evi_robustness      | evi_median_matchedcontrols.png, evi_mean_matchedcontrols.png           | Supp. Fig. 2: matched control locations |
|                     | evi_median_rand_inference.png                                          | Supp. Fig. 3: placebo-timing test |
| trees_crops         | trees_crops.tex, trees_crops_post_coefficients.png                     | Supp. Table 5 and Fig. 3: trees and crops |
| channel             | nlivestock_trees_crops.tex                                             | Supp. Table 13: change in livestock on change in trees / crops |
| health              | health.tex, health_post_coefficients.png                               | Supp. Table 6 and Fig. 4: food-security outcomes |
| add_rob             | rob_size.tex, cross_farm_char_livestock.tex                            | Supp. Tables 9 and 14: farm-size trimming, heterogeneity by household characteristics |

Steps 2-4 write `output/figures/farm_location_map.png` (Supp. Fig. 1),
`output/tables/value_scc.tex` (carbon valuation of the Methods),
`output/python_replication.csv` and `output/python_vs_stata.csv`. The `output/`
folder is shipped with the results of exactly this demo run (log, run summary,
tables, figures, comparison files), so the expected output can be inspected
before running anything. Stata's batch mode also writes a copy of the log,
`run_demo.log`, in the folder it was launched from. Step 4 prints the Python
estimates next to the Stata values; the OLS, PPML and Sun & Abraham columns
agree to the printed precision.

---------------------------------------------------------------------------------

## 4. Instructions for use (your own data)

The Stata programs read the two data files by name from `$data_dir/final/`.
To run them on other data:

1. Prepare files with the same names and variables as the demo files. The
   complete variable lists, types and the construction of the sample flags are
   documented in `code/01_make_pseudo_data.py` (section 1 for the household panel,
   section 3 for the satellite panel). The variables the estimations depend on:
   * `TREES_data_pseudo.dta` (household panel): the 88 variables listed in
     `KEEP_HH` of `code/01_make_pseudo_data.py`, i.e. every household variable the
     scripts touch: identifiers and timing (`forestgardenid`, `iso_trees`,
     `project`, `proj_area`, `projectauto`, `time_from_start`, `year_survey`,
     `month_survey`, `treat_year`, `year_cohort`); fixed-effect ids (`fe_farm`,
     `fe_proj_area`, `fe_projyr`, `fe_country`, `fe_year`); the sample flags
     `insample`, `n_obs`, `gap_ok`, `gap_ok_fl`, `first_obs`, `last_obs`,
     `post`; outcomes (`havelivestock`, `tot_animal`, `n_*count`,
     `tot_animal_ha`, `wtot_animal_ha`, `w_animals_value`, `trees`, `w_trees`,
     `w_marketableproducts`, `w_foodcrops`, `trees_per_ha`, `w1_*_ha`,
     `nofoodid`, `hfiascatid4`, `hdds_id3`, `hfias`, `hfiascategory`);
     household characteristics used in the summary and comparison tables
     (`size`, `baselinesize`, `npeople`, `householdfemales`, `n_childrennumber`,
     `n_childrenunder5`, `householdhighestschoollevel`, `farmingmainsource`,
     `childreninschool`, `householdchildren`, `worriedfood`, `fewermeals`,
     `slepthungry`, `noeating24hrs`, `guestsstaylastnight`); the local-GDP
     controls (`w5_gcpc_2017usd`, `w5ch_gcpc_2017usd`, `ch_w5_gcpc_2017usd`,
     `w1_gcpc_2017usd`, `w1ch_gcpc_2017usd`, `l1_gcpc_2017usd`, `gcp_2017usd`);
     the baseline (first-survey) controls `wb_*` and `b_n_childrennumber`; and
     the first-to-last changes `ch_tot_animal` and `wch_*`;
   * `EVI_farms_pseudo.dta` (satellite panel): `forestgardenid`, `fe_farm`,
     `year`, `calyear_treat` (0 = never treated), `never_treat`, `treat_year`,
     `treat_yearid5`-`treat_yearid13` (event time -3 to +5), `TREES_sample`,
     `projunit_calyear_treat`, `iso_trees`, `fe_ctryyr`, `fe_cell`, `fe_proj_area`,
     `w1_evi_median`, `w1_evi_mean`, `w5_evi_median`, `gcpc_2017usd_w1`,
     `gcpc_2017usd_w5`, `ch_gcpc_2017usd_w1`, `lat`, `lon` (map only).
2. Point `run_demo.do` to the data (`global data_dir`) and output folders, or set
   the three globals interactively and run `do code/03_main_estimates.do`
   followed by the program names you need (e.g. `baseline`, `evi_analysis`).
3. For the Python scripts edit `.env` (`DATA_DIR`, `OVERLEAF_FIGURES`,
   `OVERLEAF_TABLES`; relative paths are resolved from the repository root).
4. The Google Earth Engine scripts take a farm table (id, country, lat, lon,
   treatment year) uploaded as a GEE asset; the asset path, run mode
   (`control` / `treated`) and export folder are set at the top of each file.

---------------------------------------------------------------------------------

## 5. Reproduction of the manuscript's results

With the restricted data in place of the demo files (`TREES_data_pseudo.dta`
replaced by the household panel, `EVI_farms_pseudo.dta` by the satellite
panel), `run_demo.do` reproduces every table and figure of the manuscript and
Supplementary Information that is produced in Stata (the mapping in section 3);
`04_maps.py` reproduces Supplementary Fig. 1 and `05_value_scc.py` the
carbon-valuation table of the Methods (the latter needs no restricted data and
reproduces the manuscript's numbers as shipped). The upstream data-construction
scripts (aggregation of the GEE exports, survey cleaning, merging of the
local-GDP grid, construction of the sample flags) run only on the restricted raw
files and are available from the corresponding author on request, together with
the GEE asset identifiers.

Known details worth stating for anyone re-running the code:

* `esttab` tables report **t statistics** in parentheses (esttab's default).
* The placebo-timing test (`evi_rob_randinf`) is shipped with 5 placebo draws
  (`local R = 5` in `03_main_estimates.do`) so that the whole demo runs in
  under two minutes; the manuscript uses 1,000 draws. Set `local R = 1000` to
  reproduce Supp. Fig. 3 (about 40 minutes on a laptop with the demo data).

---------------------------------------------------------------------------------

## License

Code: MIT License (see `LICENSE`). Synthetic demo data: CC0 1.0. Natural Earth
boundaries: public domain. When using the code please cite the manuscript and
the repository (https://github.com/claudiorizzi/natural_capital_restoration).
