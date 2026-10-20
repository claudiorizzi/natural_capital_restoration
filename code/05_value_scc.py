# -*- coding: utf-8 -*-
"""
TREES SCC and Carbon Credit Valuation Script

What this script does
---------------------
1. Computes the present value (PV) of climate benefits for all TREES projects
   under multiple Social Cost of Carbon (SCC) values and survival rates.

2. Separately computes the PV of carbon-credit revenues for a single farmer
   completing the program, using a carbon price of $7.70 per tonne.

Key fixes relative to the original script
-----------------------------------------
- Handles multiple SCC values and survival assumptions automatically
- Fixes the 10-year horizon issue for the single-farmer calculation
- Uses a constant carbon price for the farmer carbon-credit estimate
- Prints all requested estimates without manual editing

Changes in v4
-------------
- Sequestration semantics switched to "trees alive end-of-year": a cohort
  planted in year tau contributes P_tau * survival^(t - tau + 1) * phi in
  year t (so newly planted trees do not get a full year of phi).
- CO2 totals are split into historical (year < base_year) and future
  (base_year <= year <= co2_end_year) blocks rather than a single hardcoded
  "co2_to_2035" column.
- All-TREES PV results are reported in USD millions for readability.

Changes in v5
-------------
- All-TREES valuation horizon shortened from 2100 to 2035 so the PV
  matches the manuscript's stated t_max = 2035 (conservative horizon).
  Previously PV summed benefits through 2100 while only the CO2 tonnage
  columns stopped at 2035.
- (2026-09-17 update) Single-farmer carbon-credit revenues (Estimates 2 and
  2b) are now computed at two discount rates, replacing the single 12.6%:
  10% (lower bound on PV, following Busch et al. 2019) and 2% (upper bound,
  the social discount rate used in SCC estimates). Panel B of value_scc.tex
  gains a "Discount Rate" column; the compact summary reports the rate too.
  The all-TREES SCC valuation is unchanged at 2%.
- (2026-09-17 update) Annual survival rates changed from 0.926 / 0.9756 to
  0.933 / 0.978, which imply the 50% / 80% ten-year survival stated in the
  manuscript (the old values implied 46% / 78%).

"""

# Replication package: identical to the authors' Code/7_value_scc_v5.py except that the
# paths read from .env are resolved relative to the package root (so the
# shipped .env with relative paths works from any working directory).

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
from dotenv import load_dotenv


# =========================
# Configuration
# =========================

@dataclass(frozen=True)
class SCCConfig:
    use_time_varying: bool = True
    scc_base: float = 185.0   # USD / tCO2 at base_year
    scc_growth: float = 0.02  # annual real SCC growth if time-varying


@dataclass(frozen=True)
class TreesConfig:
    start_year: int = 1989
    base_year: int = 2025
    proj_years: int = 10
    end_year: int = 2035

    survival: float = 0.933      # annual survival probability (50% over 10 years)
    phi: float = 0.02            # tCO2 per surviving tree per year
    planting_growth: float = 0.1543
    discount_rate: float = 0.02
    future_only: bool = False


# =========================
# Core functions
# =========================

def build_planting_path(
    cfg: TreesConfig,
    milestones: Dict[int, float],
    interpolation: str = "linear",
) -> pd.DataFrame:
    """
    Builds annual planting counts P_t from milestone cumulative totals,
    then projects forward planting using cfg.planting_growth for cfg.proj_years.

    Returns a DataFrame with columns:
      - year
      - P_t   (trees planted in year)
      - C_cum (cumulative planted up to year)
    """
    years_hist = np.arange(cfg.start_year, cfg.base_year + 1)
    df_hist = pd.DataFrame({"year": years_hist, "C_cum": np.nan})

    for y, cum in milestones.items():
        if cfg.start_year <= y <= cfg.base_year:
            df_hist.loc[df_hist["year"] == y, "C_cum"] = float(cum)

    df_hist = df_hist.sort_values("year").reset_index(drop=True)

    # Ensure start year exists
    if pd.isna(df_hist.loc[df_hist["year"] == cfg.start_year, "C_cum"]).all():
        df_hist.loc[df_hist["year"] == cfg.start_year, "C_cum"] = 0.0

    # Interpolate between milestone years
    df_hist["C_cum"] = df_hist["C_cum"].interpolate(method=interpolation)

    # Forward/backward fill edge cases if needed
    df_hist["C_cum"] = df_hist["C_cum"].ffill().bfill()

    # Planting flow = annual change in cumulative planted
    df_hist["P_t"] = df_hist["C_cum"].diff()
    df_hist.loc[df_hist["year"] == cfg.start_year, "P_t"] = df_hist.loc[
        df_hist["year"] == cfg.start_year, "C_cum"
    ].iloc[0]

    # Future planting from base_year + 1 through base_year + proj_years
    last_hist_P = float(df_hist.loc[df_hist["year"] == cfg.base_year, "P_t"].iloc[0])
    years_proj = np.arange(cfg.base_year + 1, cfg.base_year + cfg.proj_years + 1)

    P_proj: List[float] = []
    P_prev = last_hist_P
    for _y in years_proj:
        P_new = P_prev * (1.0 + cfg.planting_growth)
        P_proj.append(P_new)
        P_prev = P_new

    df_proj = pd.DataFrame({"year": years_proj, "P_t": P_proj})
    last_cum = float(df_hist["C_cum"].iloc[-1])
    df_proj["C_cum"] = last_cum + df_proj["P_t"].cumsum()

    df = pd.concat(
        [df_hist[["year", "P_t", "C_cum"]], df_proj[["year", "P_t", "C_cum"]]],
        ignore_index=True,
    ).sort_values("year").reset_index(drop=True)

    return df


def expand_year_grid(cfg: TreesConfig, planting_df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands to a full year grid [start_year, end_year] and merges planting flows.
    Years beyond the last planting year get P_t = 0.
    """
    years_full = np.arange(cfg.start_year, cfg.end_year + 1)
    df_full = pd.DataFrame({"year": years_full})
    df_full = df_full.merge(planting_df[["year", "P_t"]], on="year", how="left")
    df_full["P_t"] = df_full["P_t"].fillna(0.0)
    return df_full


def compute_annual_sequestration(cfg: TreesConfig, df_full: pd.DataFrame) -> pd.DataFrame:
    """
    Computes annual sequestration with "trees alive end-of-year" semantics:
      Q_t = sum_{tau <= t} P_tau * survival^(t - tau + 1) * phi
    A cohort planted in year tau contributes survival * phi in its planting
    year (not a full phi), reflecting one year of mortality before counting.
    """
    years = df_full["year"].to_numpy()
    P = df_full["P_t"].to_numpy()

    Q_t = np.zeros(len(df_full), dtype=float)

    for i, t in enumerate(years):
        ages = t - years[: i + 1] + 1
        survivors = P[: i + 1] * (cfg.survival ** ages)
        Q_t[i] = survivors.sum() * cfg.phi

    out = df_full.copy()
    out["Q_t"] = Q_t
    return out


def compute_scc_path(cfg: TreesConfig, scc_cfg: SCCConfig, df_full: pd.DataFrame) -> pd.DataFrame:
    """
    Adds SCC_t column.
    """
    out = df_full.copy()
    if scc_cfg.use_time_varying:
        out["SCC_t"] = scc_cfg.scc_base * (1.0 + scc_cfg.scc_growth) ** (
            out["year"] - cfg.base_year
        )
    else:
        out["SCC_t"] = scc_cfg.scc_base
    return out


def compute_pv(cfg: TreesConfig, df_full: pd.DataFrame) -> Tuple[float, Optional[float], float, pd.DataFrame]:
    """
    Computes annual benefits B_t = Q_t * SCC_t and discounted PV to base_year.
    """
    out = df_full.copy()
    out["B_t"] = out["Q_t"] * out["SCC_t"]
    out["disc_factor"] = (1.0 + cfg.discount_rate) ** (out["year"] - cfg.base_year)
    out["B_t_PV"] = out["B_t"] / out["disc_factor"]

    if cfg.future_only:
        mask_future = out["year"] >= cfg.base_year
        pv_future = float(out.loc[mask_future, "B_t_PV"].sum())
        return pv_future, None, pv_future, out

    mask_hist = out["year"] < cfg.base_year
    mask_future = out["year"] >= cfg.base_year

    pv_hist = float(out.loc[mask_hist, "B_t_PV"].sum())
    pv_future = float(out.loc[mask_future, "B_t_PV"].sum())
    pv_total = float(out["B_t_PV"].sum())

    return pv_total, pv_hist, pv_future, out


def run_valuation(
    cfg: TreesConfig,
    scc_cfg: SCCConfig,
    milestones: Dict[int, float],
) -> Tuple[float, Optional[float], float, pd.DataFrame]:
    """
    End-to-end runner.
    """
    planting = build_planting_path(cfg, milestones)
    df_full = expand_year_grid(cfg, planting)
    df_full = compute_annual_sequestration(cfg, df_full)
    df_full = compute_scc_path(cfg, scc_cfg, df_full)
    return compute_pv(cfg, df_full)


# =========================
# Reporting helpers
# =========================

def summarize_scenario(
    label: str,
    cfg: TreesConfig,
    scc_cfg: SCCConfig,
    milestones: Dict[int, float],
    co2_end_year: int = 2035,
) -> pd.Series:
    """
    Runs one scenario and returns a tidy summary row.

    CO2 totals are split into:
      - co2_hist:   cumulative Q_t for years strictly before base_year
      - co2_future: cumulative Q_t for base_year through co2_end_year (inclusive)
    """
    pv_total, pv_hist, pv_future, df_vals = run_valuation(cfg, scc_cfg, milestones)

    mask_hist = df_vals["year"] < cfg.base_year
    mask_future = (df_vals["year"] >= cfg.base_year) & (df_vals["year"] <= co2_end_year)
    co2_hist_tons = float(df_vals.loc[mask_hist, "Q_t"].sum())
    co2_future_tons = float(df_vals.loc[mask_future, "Q_t"].sum())

    row = {
        "scenario": label,
        "survival": cfg.survival,
        "scc_base": scc_cfg.scc_base,
        "time_varying_scc": scc_cfg.use_time_varying,
        "scc_growth": scc_cfg.scc_growth,
        "discount_rate": cfg.discount_rate,
        "future_only": cfg.future_only,
        "co2_end_year": co2_end_year,
        "pv_historical_usd": pv_hist,
        "pv_future_usd": pv_future,
        "pv_total_usd": pv_total,
        "co2_hist_tons": co2_hist_tons,
        "co2_future_tons": co2_future_tons,
        "co2_hist_mt": co2_hist_tons / 1e6,
        "co2_future_mt": co2_future_tons / 1e6,
    }
    return pd.Series(row)


# =========================
# LaTeX output
# =========================

def write_latex_tables(
    df_all: pd.DataFrame,
    df_farmer: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Writes a combined LaTeX table with two panels to output_path:
      Panel A: All TREES projects (PV in USD millions)
      Panel B: Single farmer carbon credits (PV in USD)
    Booktabs style — requires \\usepackage{booktabs} in the LaTeX preamble.
    """
    lines: List[str] = []
    lines.append(r"\begin{table}[!htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Present Value of Climate Benefits from TREES}")
    lines.append(r"\label{tab:value_scc}")
    lines.append(r"\small")
    lines.append("")

    # Panel A: All TREES
    lines.append(r"\textbf{Panel A: All TREES Projects (1989--2035)}")
    lines.append("")
    lines.append(r"\begin{tabular}{rcrrrrr}")
    lines.append(r"\toprule")
    lines.append(
        r"SCC & Survival & PV Hist & PV Future & PV Total & CO$_2$ Hist & CO$_2$ Future \\"
    )
    lines.append(
        r"(\$/tCO$_2$) & & (USD M) & (USD M) & (USD M) & (Mt) & (Mt) \\"
    )
    lines.append(r"\midrule")
    for _, row in df_all.iterrows():
        cells = [
            f"\\${row['scc_base']:,.0f}",
            f"{row['survival']:.4f}",
            f"{row['pv_historical_musd']:,.2f}",
            f"{row['pv_future_musd']:,.2f}",
            f"{row['pv_total_musd']:,.2f}",
            f"{row['co2_hist_mt']:,.2f}",
            f"{row['co2_future_mt']:,.2f}",
        ]
        lines.append(" & ".join(cells) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append("")
    lines.append(r"\vspace{1em}")
    lines.append("")

    # Panel B: Single farmer
    lines.append(
        r"\textbf{Panel B: Single Farmer Carbon Credits (2025--2034, price = \$7.70/tCO$_2$)}"
    )
    lines.append("")
    lines.append(r"\begin{tabular}{ccrr}")
    lines.append(r"\toprule")
    lines.append(r"Discount Rate & Survival & PV Total & CO$_2$ Future \\")
    lines.append(r" & & (USD) & (tonnes) \\")
    lines.append(r"\midrule")
    for _, row in df_farmer.iterrows():
        cells = [
            f"{100 * row['discount_rate']:.0f}\\%",
            f"{row['survival']:.4f}",
            f"{row['pv_total_usd']:,.2f}",
            f"{row['co2_future_tons']:,.2f}",
        ]
        lines.append(" & ".join(cells) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append("")
    lines.append(r"\end{table}")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# =========================
# Main
# =========================

def main() -> None:
    # -------------------------------------------------
    # Estimate 1: All TREES projects since founding
    # -------------------------------------------------
    milestones_all_trees = {
        1989: 0,
        2003: 30_000_000,
        2015: 100_000_000,
        2025: 420_000_000,
    }

    # Annual survival consistent with 50% / 80% ten-year survival
    # (0.933^10 = 0.50, 0.978^10 = 0.80), matching the manuscript
    survival_rates = [0.933, 0.978]
    scc_values = [185.0, 77.0]

    results_all = []

    for scc in scc_values:
        for survival in survival_rates:
            cfg = TreesConfig(
                start_year=1989,
                base_year=2025,
                proj_years=10,
                end_year=2035,   # t_max = 2035, matching the manuscript's conservative horizon
                survival=survival,
                phi=0.02,
                planting_growth=0.1543,
                discount_rate=0.02,
                future_only=False,
            )

            scc_cfg = SCCConfig(
                use_time_varying=True,
                scc_base=scc,
                scc_growth=0.02,
            )

            label = f"All TREES | SCC=${scc:,.0f} | survival={survival}"
            results_all.append(
                summarize_scenario(label, cfg, scc_cfg, milestones_all_trees, co2_end_year=2035)
            )

    df_results_all = pd.DataFrame(results_all)

    # Convert PV columns from USD to USD millions for readability
    df_results_all["pv_historical_musd"] = df_results_all["pv_historical_usd"] / 1e6
    df_results_all["pv_future_musd"] = df_results_all["pv_future_usd"] / 1e6
    df_results_all["pv_total_musd"] = df_results_all["pv_total_usd"] / 1e6

    print("\n" + "=" * 80)
    print("ALL TREES PROJECTS: PRESENT VALUE OF CLIMATE BENEFITS (USD millions)")
    print("=" * 80)
    cols_all = [
        "scenario",
        "pv_historical_musd",
        "pv_future_musd",
        "pv_total_musd",
        "co2_hist_mt",
        "co2_future_mt",
    ]
    print(
        df_results_all[cols_all].to_string(
            index=False, float_format=lambda x: f"{x:,.2f}"
        )
    )

    # -------------------------------------------------
    # Estimate 2: Single farmer carbon-credit revenues
    # -------------------------------------------------
    # Average of number of trees 3 years from beginning of project: 1,069
    milestones_farmer = {
        2025: 1069,
    }

    # Discount-rate bounds for the farmer-level offset revenues (farmer level
    # only; the all-TREES SCC valuation above keeps the 2% social rate):
    #   lower bound on PV: 10%, following Busch et al. (2019)
    #   upper bound on PV:  2%, the social discount rate used in SCC estimates
    farmer_discount_rates = [0.10, 0.02]

    results_farmer = []

    for rate in farmer_discount_rates:
        for survival in survival_rates:
            cfg_farmer = TreesConfig(
                start_year=2025,
                base_year=2025,
                proj_years=0,          # no new planting after initial stock
                end_year=2034,         # 10-year horizon: 2025-2034 inclusive
                survival=survival,
                phi=0.02,
                planting_growth=0.0,
                discount_rate=rate,
                future_only=True,
            )

            carbon_price_cfg = SCCConfig(
                use_time_varying=False,  # fixed carbon price of $7.70
                scc_base=7.70,
                scc_growth=0.0,
            )

            label = (f"Single farmer carbon credits | price=$7.70 | "
                     f"r={rate:.0%} | survival={survival}")
            results_farmer.append(
                summarize_scenario(label, cfg_farmer, carbon_price_cfg, milestones_farmer, co2_end_year=2034)
            )

    df_results_farmer = pd.DataFrame(results_farmer)

    print("\n" + "=" * 80)
    print("SINGLE FARMER: PRESENT VALUE OF CARBON-CREDIT REVENUES (USD)")
    print("=" * 80)
    cols_farmer = [
        "scenario",
        "pv_total_usd",
        "co2_future_tons",
    ]
    print(
        df_results_farmer[cols_farmer].to_string(
            index=False, float_format=lambda x: f"{x:,.2f}"
        )
    )

    # -------------------------------------------------
    # Estimate 2b: Single farmer carbon-credit revenues WITH TREE REPLACEMENT
    # -------------------------------------------------
    # Assumption: the farmer replaces every tree that dies, so the stock of
    # trees doing carbon sequestration stays constant at the initial 1,069
    # over the whole 2025-2034 horizon. Mechanically this is survival = 1.0
    # in the cohort model (each cohort contributes a full phi every year).
    # The cost of replacement seedlings is NOT netted out of the revenue.
    # Same 10% / 2% discount-rate bounds as Estimate 2.
    results_farmer_repl = []

    for rate in farmer_discount_rates:
        cfg_farmer_repl = TreesConfig(
            start_year=2025,
            base_year=2025,
            proj_years=0,          # no new planting beyond the (replaced) stock
            end_year=2034,         # 10-year horizon: 2025-2034 inclusive
            survival=1.0,          # replacement keeps the stock constant
            phi=0.02,
            planting_growth=0.0,
            discount_rate=rate,
            future_only=True,
        )

        carbon_price_cfg_repl = SCCConfig(
            use_time_varying=False,  # fixed carbon price of $7.70
            scc_base=7.70,
            scc_growth=0.0,
        )

        label_repl = (f"Single farmer carbon credits | price=$7.70 | "
                      f"r={rate:.0%} | replacement (constant stock)")
        results_farmer_repl.append(
            summarize_scenario(
                label_repl, cfg_farmer_repl, carbon_price_cfg_repl,
                milestones_farmer, co2_end_year=2034,
            )
        )

    df_results_farmer_repl = pd.DataFrame(results_farmer_repl)

    print("\n" + "=" * 80)
    print("SINGLE FARMER WITH REPLACEMENT (CONSTANT 1,069-TREE STOCK): PV OF CARBON-CREDIT REVENUES (USD)")
    print("=" * 80)
    print(
        df_results_farmer_repl[cols_farmer].to_string(
            index=False, float_format=lambda x: f"{x:,.2f}"
        )
    )

    # Include the replacement scenario in the LaTeX Panel B and the compact
    # summary below (its survival column shows 1.0000 = constant stock).
    df_results_farmer = pd.concat(
        [df_results_farmer, df_results_farmer_repl], ignore_index=True
    )

    # Compact table for easy copy/paste. PV is reported per-row in its native
    # unit (USD millions for All TREES, USD for the single farmer) with the
    # unit shown in its own column to avoid mixing scales.
    print("\n" + "=" * 80)
    print("COMPACT SUMMARY TABLE")
    print("=" * 80)

    compact_rows = []

    for _, row in df_results_all.iterrows():
        compact_rows.append({
            "type": "All TREES",
            "SCC_or_price": row["scc_base"],
            "discount_rate": row["discount_rate"],
            "survival": row["survival"],
            "PV_total": row["pv_total_usd"] / 1e6,
            "units": "USD millions",
        })

    for _, row in df_results_farmer.iterrows():
        compact_rows.append({
            "type": "Single farmer credits",
            "SCC_or_price": 7.70,
            "discount_rate": row["discount_rate"],
            "survival": row["survival"],
            "PV_total": row["pv_total_usd"],
            "units": "USD",
        })

    df_compact = pd.DataFrame(compact_rows)
    print(
        df_compact.to_string(
            index=False, float_format=lambda x: f"{x:,.2f}"
        )
    )

    # -------------------------------------------------
    # Write LaTeX table to Overleaf Tables folder
    # -------------------------------------------------
    project_root = Path(__file__).resolve().parents[1]
    load_dotenv(project_root / ".env")
    overleaf_tables = os.environ.get("OVERLEAF_TABLES")
    if not overleaf_tables:
        raise RuntimeError(
            "OVERLEAF_TABLES is not set in the project's .env file."
        )
    output_path = (project_root / overleaf_tables).resolve() / "value_scc.tex"
    output_path.parent.mkdir(parents=True, exist_ok=True)   # replication package: create output/tables if absent
    write_latex_tables(df_results_all, df_results_farmer, output_path)
    print(f"\nWrote LaTeX table to {output_path}")


if __name__ == "__main__":
    main()
