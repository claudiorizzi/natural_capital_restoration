# -*- coding: utf-8 -*-
"""
06_python_replication.py - independent Python re-implementation of the paper's core
estimators, run on the demo (pseudo) data, used to cross-check the Stata output.

Why this exists
---------------
Nature's code policy asks that reviewers be able to run the software. The main
analysis is written in Stata with user-written packages (reghdfe, ppmlhdfe,
jwdid, eventstudyinteract). This script re-implements the estimators that
produce the paper's headline tables with numpy/pandas only, so that
  (i)  the demo data can be analysed without a Stata licence, and
  (ii) the Stata pipeline can be verified against an independent implementation
       (see `compare_with_stata()`; the comparison is printed when the Stata
       demo has been run and output/tables/*.tex exist).

What is replicated (programs in 03_main_estimates.do)
----------------------------------------------------------
A. `baseline`  -> Table own_n_livestock (4 columns)
     (1)-(2) OLS   havelivestock ~ post [+ GCP]      FE: month_survey, fe_farm; cluster fe_proj_area
     (3)-(4) PPML  tot_animal    ~ post [+ GCP]      FE: month_survey, fe_farm; cluster fe_proj_area
B. `trees_crops` -> Table trees_crops (6 PPML columns, same design)
C. `health`      -> Table health (3 OLS columns, same design, + post2019 in cols 1-2)
D. `evi_analysis` -> Table evi_median_mean and Fig. evi_*_combined
     event-study ATTs for median and mean EVI (all clustered by project area, fe_proj_area) with
       - two-way FE event study (the `reghdfe` sample-definition helper),
       - Wooldridge (2021) extended TWFE ("stacked DiD", jwdid in Stata),
       - Sun & Abraham (2021) interaction-weighted estimator (eventstudyinteract).
   Every sample restriction of the Stata code is reproduced: the [-3,+4] event
   window and the reference period t = -1.

Estimator details (match reghdfe / ppmlhdfe conventions)
--------------------------------------------------------
* Fixed effects are removed by alternating projections (within transformation),
  weighted by the IRLS weights for PPML.
* Standard errors are cluster-robust with reghdfe's small-sample adjustment
  q = (N-1)/(N-K) * G/(G-1), where K counts the regressors plus the absorbed
  fixed-effect levels that are NOT nested inside the cluster variable (minus the
  redundant levels, i.e. the number of mobility groups for the second FE).
* PPML follows ppmlhdfe: IRLS with HDFE, observations in fixed-effect groups whose
  outcome is always zero are dropped (they are perfectly predicted), then singleton
  groups are dropped iteratively.
* The ETWFE implementation enters the GCP covariate linearly (jwdid by default
  also interacts covariates with cohort and period dummies), and the Sun & Abraham
  variance omits the (tiny) sampling-variance term of the cohort shares that
  eventstudyinteract adds. These are the only known deviations; both are stated
  in the comparison output.

Usage (from the repository root):
    uv run python code/06_python_replication.py
Output:
    output/python_replication.csv     all estimates (spec, outcome, term, b, se, N)
    console                           tables + comparison with Stata .tex files
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from dotenv import load_dotenv

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
load_dotenv(ROOT / ".env")
DATA = (ROOT / os.environ.get("DATA_DIR", "data")).resolve()
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)
TABLES = (ROOT / os.environ.get("OVERLEAF_TABLES", "output/tables")).resolve()

# =======================================================
# 1. Fixed-effect regression toolkit
# =======================================================

def _codes(a):
    """integer group codes 0..L-1 for an array-like"""
    return pd.factorize(pd.Series(np.asarray(a)), sort=True)[0]


def demean(Y, fes, w=None, tol=1e-9, maxiter=5000):
    """Remove fixed effects from the columns of Y (n x p) by alternating
    projections (Gauss-Seidel over the FE groups), optionally weighted."""
    Y = np.array(Y, dtype=float, copy=True)
    if Y.ndim == 1:
        Y = Y[:, None]
    w = np.ones(len(Y)) if w is None else np.asarray(w, float)
    fes = [np.asarray(f) for f in fes]
    sizes = [np.bincount(f, weights=w) for f in fes]
    for _ in range(maxiter):
        delta = 0.0
        for f, s in zip(fes, sizes):
            for j in range(Y.shape[1]):
                m = np.bincount(f, weights=w * Y[:, j]) / s
                adj = m[f]
                Y[:, j] -= adj
                delta = max(delta, np.abs(adj).max())
        if delta < tol:
            break
    return Y


def mobility_groups(f1, f2):
    """number of connected components of the bipartite graph (f1 level, f2 level)"""
    n1 = f1.max() + 1
    parent = np.arange(n1 + f2.max() + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in set(zip(f1.tolist(), (f2 + n1).tolist())):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(x) for x in range(len(parent))})


def df_absorbed(fes, cluster):
    """reghdfe's e(df_a): total FE levels, minus the redundant levels (one
    mobility group per additional FE), minus every level of an FE that is nested
    inside the cluster variable (verified against e(df_a_initial),
    e(df_a_redundant), e(df_a_nested) on the demo data)."""
    df = 0
    for i, f in enumerate(fes):
        levels = f.max() + 1
        df += levels - (mobility_groups(fes[0], f) if i > 0 else 0)
        nested = pd.Series(cluster).groupby(f).nunique().max() == 1
        if nested:
            df -= levels
    return df


def cluster_vcov(Xd, score, bread, cluster, K, dof_adjust=True):
    """q * bread (sum_g s_g s_g') bread. reghdfe: q = (N-1)/(N-K) * G/(G-1);
    ppmlhdfe (as Stata's poisson): q = G/(G-1) only (dof_adjust=False)."""
    n = len(score)
    G = len(np.unique(cluster))
    sg = pd.DataFrame(Xd * score[:, None]).groupby(cluster).sum().to_numpy()
    S = sg.T @ sg
    q = ((n - 1) / (n - K) if dof_adjust else 1.0) * G / (G - 1)
    return q * bread @ S @ bread, G


def drop_singletons(df, fe_cols):
    """iteratively drop observations that are alone in any FE group (reghdfe default)"""
    while True:
        n0 = len(df)
        for c in fe_cols:
            df = df[df.groupby(c)[c].transform("size") > 1]
        if len(df) == n0:
            return df


def ols_hdfe(df, y, xs, fe_cols, cluster):
    """reghdfe y xs, absorb(fe_cols) vce(cluster cluster)"""
    d = df.dropna(subset=[y] + xs + fe_cols + [cluster]).copy()
    d = drop_singletons(d, fe_cols)
    fes = [_codes(d[c]) for c in fe_cols]
    cl = _codes(d[cluster])
    Z = demean(d[[y] + xs].to_numpy(float), fes)
    yd, Xd = Z[:, 0], Z[:, 1:]
    XtX = Xd.T @ Xd
    bread = np.linalg.inv(XtX)
    b = bread @ Xd.T @ yd
    e = yd - Xd @ b
    K = len(xs) + df_absorbed(fes, cl)
    V, G = cluster_vcov(Xd, e, bread, cl, K)
    r2 = 1 - (e @ e) / ((d[y] - d[y].mean()) ** 2).sum()
    return dict(b=b, se=np.sqrt(np.diag(V)), V=V, N=len(d), G=G, names=xs, r2=r2, K=K)


def ppml_hdfe(df, y, xs, fe_cols, cluster, tol=1e-8, maxiter=100):
    """ppmlhdfe y xs, absorb(fe_cols) vce(cluster cluster)"""
    d = df.dropna(subset=[y] + xs + fe_cols + [cluster]).copy()
    # drop FE groups in which y is always zero (perfectly predicted), then singletons
    while True:
        n0 = len(d)
        for c in fe_cols:
            d = d[d.groupby(c)[y].transform("sum") > 0]
        d = drop_singletons(d, fe_cols)
        if len(d) == n0:
            break
    fes = [_codes(d[c]) for c in fe_cols]
    cl = _codes(d[cluster])
    yv = d[y].to_numpy(float)
    X = d[xs].to_numpy(float)
    mu = (yv + yv.mean()) / 2                      # starting values
    eta = np.log(mu)
    b = np.zeros(len(xs))
    dev_old = np.inf
    for it in range(maxiter):
        z = eta + (yv - mu) / mu                   # working dependent variable
        Z = demean(np.column_stack([z, X]), fes, w=mu)
        zd, Xd = Z[:, 0], Z[:, 1:]
        W = mu
        XtWX = Xd.T @ (Xd * W[:, None])
        b = np.linalg.solve(XtWX, Xd.T @ (W * zd))
        r = zd - Xd @ b                            # weighted-demeaned residual
        eta = np.clip(z - r, -30, 30)              # = X b + fixed effects
        mu = np.exp(eta)
        with np.errstate(divide="ignore", invalid="ignore"):
            dev = 2 * np.sum(np.where(yv > 0, yv * np.log(yv / mu), 0) - (yv - mu))
        if abs(dev - dev_old) < tol * (1 + abs(dev)):
            break
        dev_old = dev
    # sandwich on the final weighted within-transformed regressors
    Z = demean(np.column_stack([np.zeros(len(yv)), X]), fes, w=mu)
    Xd = Z[:, 1:]
    bread = np.linalg.inv(Xd.T @ (Xd * mu[:, None]))
    K = len(xs) + df_absorbed(fes, cl)
    V, G = cluster_vcov(Xd, yv - mu, bread, cl, K, dof_adjust=False)
    return dict(b=b, se=np.sqrt(np.diag(V)), V=V, N=len(d), G=G, names=xs, r2=np.nan, K=K)


# =======================================================
# 2. Data preparation mirroring the Stata programs
# =======================================================
def load_household():
    hh = pd.read_stata(DATA / "final/TREES_data_pseudo.dta", convert_categoricals=False)
    hh["w5_gcpc_2017usd"] = hh.w5_gcpc_2017usd / 1000            # GCP in $1,000 as in `baseline`
    est = hh[(hh.insample == 1) & (hh.n_obs > 1) & (hh.gap_ok == 1) & (hh.gap_ok_fl == 1)].copy()
    est["post2019"] = (est.year_survey > 2019).astype(float)
    return est


def load_evi(window=(-3, 4)):
    ev = pd.read_stata(DATA / "final/EVI_farms_pseudo.dta", convert_categoricals=False)
    ev["event_window"] = np.where(ev.calyear_treat != 0, ev.year - ev.calyear_treat, np.nan)
    ev = ev[ev.event_window.between(*window) | (ev.calyear_treat == 0)].copy()
    ev["fe_cellyr"] = ev.groupby(["fe_cell", "year"]).ngroup() + 1
    return ev


# =======================================================
# 3. Event-study estimators (EVI)
# =======================================================
EVENT_IDS = {5: -3, 6: -2, 8: 0, 9: 1, 10: 2, 11: 3, 12: 4}     # treat_yearidN -> event time N-8 ($treat_years)


def twfe_event_study(ev, y, control, fe_cols, cluster):
    """reghdfe y $treat_years control, absorb(fe_cols) vce(cluster cluster)"""
    xs = [f"treat_yearid{k}" for k in EVENT_IDS] + [control]
    res = ols_hdfe(ev, y, xs, fe_cols, cluster)
    rows = [dict(event=EVENT_IDS[k], b=res["b"][i], se=res["se"][i]) for i, k in enumerate(EVENT_IDS)]
    return pd.DataFrame(rows), res


def _cohort_period_design(ev, kind):
    """indicator columns for cohort x (calendar year | relative period) cells of
    treated farms; the reference period t = g-1 (relative -1) is omitted."""
    treated = ev.calyear_treat > 0
    cols, names, meta = [], [], []
    for g in sorted(ev.loc[treated, "calyear_treat"].unique()):
        for rel in range(-3, 5):
            if rel == -1:
                continue
            m = (treated & (ev.calyear_treat == g) & (ev.event_window == rel)).to_numpy()
            if m.sum() == 0:
                continue
            cols.append(m.astype(float)); names.append(f"g{int(g)}_e{rel}"); meta.append((int(g), rel, int(m.sum())))
    return np.column_stack(cols), names, meta


def etwfe_event_study(ev, y, control, fe_cols, cluster):
    """Wooldridge (2021) extended TWFE with never-treated controls (jwdid ..., never):
    cohort x period cells of treated farms, FE, covariate; `estat event` aggregates
    the cell ATTs by relative period with weights = number of treated farm-years."""
    D, names, meta = _cohort_period_design(ev, "etwfe")
    d = ev.copy()
    for j, nm in enumerate(names):
        d[nm] = D[:, j]
    res = ols_hdfe(d, y, names + [control], fe_cols, cluster)
    return _aggregate_by_event(res, names, meta), res


def sunab_event_study(ev, y, fe_cols, cluster):
    """Sun & Abraham (2021) interaction-weighted estimator (eventstudyinteract ...,
    cohort(calyear_treat) control_cohort(never_treat)): cohort x relative-period
    dummies for the 7 periods in $treat_years, weights = cohort shares at each
    relative period among treated observations."""
    treated = ev.calyear_treat > 0
    cols, names, meta = [], [], []
    for g in sorted(ev.loc[treated, "calyear_treat"].unique()):
        for k, rel in EVENT_IDS.items():
            m = (treated & (ev.calyear_treat == g) & (ev[f"treat_yearid{k}"] == 1)).to_numpy()
            if m.sum() == 0:
                continue
            cols.append(m.astype(float)); names.append(f"g{int(g)}_e{rel}"); meta.append((int(g), rel, int(m.sum())))
    d = ev.copy()
    for j, nm in enumerate(names):
        d[nm] = cols[j]
    res = ols_hdfe(d, y, names, fe_cols, cluster)
    return _aggregate_by_event(res, names, meta), res


def _aggregate_by_event(res, names, meta):
    """weighted average of cohort-specific coefficients by relative period;
    delta-method SE from the clustered VCE"""
    rows = []
    for rel in range(-3, 5):
        idx = [i for i, (g, e, n) in enumerate(meta) if e == rel]
        if rel == -1 or not idx:
            rows.append(dict(event=rel, b=0.0 if rel == -1 else np.nan, se=0.0 if rel == -1 else np.nan))
            continue
        w = np.zeros(len(names) + (len(res["b"]) - len(names)))
        tot = sum(meta[i][2] for i in idx)
        for i in idx:
            w[i] = meta[i][2] / tot
        rows.append(dict(event=rel, b=float(w @ res["b"]), se=float(np.sqrt(w @ res["V"] @ w))))
    return pd.DataFrame(rows)


# =======================================================
# 4. Run everything
# =======================================================
def stars(b, se):
    z = abs(b / se) if se > 0 else 0
    return "***" if z > 2.576 else "**" if z > 1.96 else "*" if z > 1.645 else ""


def fmt(b, se):
    return f"{b:8.4f}{stars(b, se):3s} ({se:.5f})"


def main():
    results = []
    est = load_household()
    print(f"household estimation sample: {len(est):,} farm-surveys, {est.fe_farm.nunique():,} farms, "
          f"{est.fe_proj_area.nunique()} project-area clusters\n")

    # ---- A. own_n_livestock ------------------------------------------------------
    print("A. Table own_n_livestock (program `baseline`)")
    specs = [("OLS", "havelivestock", ["post"], "fe_proj_area"),
             ("OLS", "havelivestock", ["post", "w5_gcpc_2017usd"], "fe_proj_area"),
             ("PPML", "tot_animal", ["post"], "fe_proj_area"),
             ("PPML", "tot_animal", ["post", "w5_gcpc_2017usd"], "fe_proj_area")]
    for i, (kind, y, xs, cl) in enumerate(specs, 1):
        f = ols_hdfe if kind == "OLS" else ppml_hdfe
        r = f(est, y, xs, ["month_survey", "fe_farm"], cl)
        print(f"  ({i}) {kind:4s} {y:14s} N={r['N']:5d} G={r['G']:3d} K={r['K']:4d}  "
              + "  ".join(f"{nm}={fmt(b, s)}" for nm, b, s in zip(r["names"], r["b"], r["se"])))
        for nm, b, s in zip(r["names"], r["b"], r["se"]):
            results.append(dict(table="own_n_livestock", col=i, estimator=kind, outcome=y, term=nm, b=b, se=s, N=r["N"], G=r["G"]))

    # ---- B. trees_crops ----------------------------------------------------------------
    print("\nB. Table trees_crops (program `trees_crops`, PPML, cluster fe_proj_area)")
    for i, y in enumerate(["w_trees", "w_marketableproducts", "w_foodcrops", "w1_trees_per_ha", "w1_mark_prod_ha", "w1_foodcrop_ha"], 1):
        r = ppml_hdfe(est, y, ["post", "w5_gcpc_2017usd"], ["month_survey", "fe_farm"], "fe_proj_area")
        print(f"  ({i}) {y:22s} N={r['N']:5d}  post={fmt(r['b'][0], r['se'][0])}  gcpc={fmt(r['b'][1], r['se'][1])}")
        for nm, b, s in zip(r["names"], r["b"], r["se"]):
            results.append(dict(table="trees_crops", col=i, estimator="PPML", outcome=y, term=nm, b=b, se=s, N=r["N"], G=r["G"]))

    # ---- C. health ------------------------------------------------------------------------
    print("\nC. Table health (program `health`, OLS, cluster fe_proj_area)")
    for i, (y, xs) in enumerate([("nofoodid", ["post", "w5_gcpc_2017usd", "post2019"]),
                                 ("hfiascatid4", ["post", "w5_gcpc_2017usd", "post2019"]),
                                 ("hdds_id3", ["post", "w5_gcpc_2017usd"])], 1):
        r = ols_hdfe(est, y, xs, ["month_survey", "fe_farm"], "fe_proj_area")
        print(f"  ({i}) {y:12s} N={r['N']:5d}  post={fmt(r['b'][0], r['se'][0])}  gcpc={fmt(r['b'][1], r['se'][1])}  R2={r['r2']:.3f}")
        for nm, b, s in zip(r["names"], r["b"], r["se"]):
            results.append(dict(table="health", col=i, estimator="OLS", outcome=y, term=nm, b=b, se=s, N=r["N"], G=r["G"]))

    # ---- D. EVI event studies ----------------------------------------------------------
    ev = load_evi()
    print(f"\nD. EVI event studies (program `evi_analysis`): {len(ev):,} farm-years, "
          f"{ev.forestgardenid.nunique():,} farms, {int((ev.calyear_treat > 0).sum()):,} treated farm-years")
    ev_med = ev.copy()
    # sample-definition helper (reghdfe) -> its estimation sample is used by jwdid and SA
    twfe_m, r_twfe = twfe_event_study(ev_med, "w1_evi_median", "gcpc_2017usd_w1", ["fe_farm", "fe_ctryyr"], "fe_proj_area")
    etw_m, r_etw = etwfe_event_study(ev_med, "w1_evi_median", "gcpc_2017usd_w5", ["fe_farm", "fe_ctryyr"], "fe_proj_area")
    sa_m, r_sa = sunab_event_study(ev_med, "w1_evi_median", ["fe_farm", "fe_cellyr"], "fe_proj_area")
    # mean EVI: jwdid on w1_evi_mean with cell-year FE; the Stata
    # SA "mean" column is run on w5_evi_median*100 (the script overwrites w1_evi_mean)
    ev_mean = ev.copy()
    twfe_mn, _ = twfe_event_study(ev_mean, "w1_evi_mean", "gcpc_2017usd_w5", ["fe_farm", "fe_ctryyr"], "fe_proj_area")
    etw_mn, _ = etwfe_event_study(ev_mean, "w1_evi_mean", "gcpc_2017usd_w5", ["fe_farm", "fe_cellyr"], "fe_proj_area")
    ev_mean["w1_evi_mean_asinstata"] = ev_mean.w5_evi_median * 100
    sa_mn, _ = sunab_event_study(ev_mean, "w1_evi_mean_asinstata", ["fe_farm", "fe_cellyr"], "fe_proj_area")

    tab = twfe_m.rename(columns={"b": "TWFE_b", "se": "TWFE_se"})
    for nm, dfr in [("ETWFE", etw_m), ("SunAb", sa_m), ("TWFE_mean", twfe_mn), ("ETWFE_mean", etw_mn), ("SunAb_mean", sa_mn)]:
        tab = tab.merge(dfr.rename(columns={"b": f"{nm}_b", "se": f"{nm}_se"}), on="event", how="outer")
    pd.set_option("display.width", 200)
    print(f"  median: TWFE N={r_twfe['N']:,} G={r_twfe['G']}  ETWFE N={r_etw['N']:,}  SA N={r_sa['N']:,}")
    print(tab.round(3).to_string(index=False))
    for _, row in tab.iterrows():
        for nm in ["TWFE", "ETWFE", "SunAb", "TWFE_mean", "ETWFE_mean", "SunAb_mean"]:
            results.append(dict(table="evi_median_mean", col=nm, estimator=nm.split("_")[0],
                                outcome="w1_evi_mean" if nm.endswith("mean") else "w1_evi_median",
                                term=f"event{int(row.event):+d}", b=row[f"{nm}_b"], se=row[f"{nm}_se"], N=np.nan, G=np.nan))

    out = pd.DataFrame(results)
    out.to_csv(OUT / "python_replication.csv", index=False)
    print(f"\nestimates written to {OUT / 'python_replication.csv'}")
    compare_with_stata(out)


# =======================================================
# 5. Comparison with the Stata tables (if the Stata demo has run)
# =======================================================
def parse_esttab(path):
    """return {row label: ([b...], [t...])} from an esttab booktabs table. The
    paper's esttab calls use the default parentheses content, i.e. t statistics."""
    txt = path.read_text(encoding="utf-8").splitlines()
    rows, i = {}, 0
    while i < len(txt):
        line = txt[i]
        if "&" in line and not line.lstrip().startswith(("\\", "&")) and "multicolumn" not in line:
            cells = [c.strip() for c in line.rstrip("\\").split("&")]
            label = cells[0]
            try:
                b = [float(re.sub(r"\\sym\{[^}]*\}", "", c)) if c else np.nan for c in cells[1:]]
            except ValueError:
                i += 1; continue
            se = []
            if i + 1 < len(txt) and txt[i + 1].lstrip().startswith("&"):
                se = [float(c.strip().strip("()")) if c.strip() else np.nan for c in txt[i + 1].rstrip("\\").split("&")[1:]]
                i += 1
            rows[label] = (b, se)
        i += 1
    return rows


def compare_with_stata(py):
    checks = []
    f = TABLES / "own_n_livestock.tex"
    if f.exists():
        st = parse_esttab(f)
        for label, term in [("1(Post 1st survey)", "post"), ("GCP per capita (\\textdollar 1,000)", "w5_gcpc_2017usd")]:
            if label in st:
                for col, (b, se) in enumerate(zip(*st[label]), 1):
                    q = py[(py.table == "own_n_livestock") & (py.col == col) & (py.term == term)]
                    if len(q) and not np.isnan(b):
                        checks.append(("own_n_livestock", col, term, b, se, q.b.iloc[0], q.se.iloc[0]))
    f = TABLES / "trees_crops.tex"
    if f.exists():
        st = parse_esttab(f)
        for label, term in [("1(Post 1st survey)", "post"), ("GPC per capita (\\textdollar 1,000)", "w5_gcpc_2017usd")]:
            if label in st:
                for col, (b, se) in enumerate(zip(*st[label]), 1):
                    q = py[(py.table == "trees_crops") & (py.col == col) & (py.term == term)]
                    if len(q): checks.append(("trees_crops", col, term, b, se, q.b.iloc[0], q.se.iloc[0]))
    f = TABLES / "health.tex"
    if f.exists():
        st = parse_esttab(f)
        for label, term in [("1(Post 1st survey)", "post"), ("GCP per capita (\\textdollar 1,000)", "w5_gcpc_2017usd")]:
            if label in st:
                for col, (b, se) in enumerate(zip(*st[label]), 1):
                    q = py[(py.table == "health") & (py.col == col) & (py.term == term)]
                    if len(q): checks.append(("health", col, term, b, se, q.b.iloc[0], q.se.iloc[0]))
    f = TABLES / "evi_median_mean.tex"
    if f.exists():
        st = parse_esttab(f)
        colmap = {1: "ETWFE", 2: "SunAb", 3: "ETWFE_mean", 4: "SunAb_mean"}
        for k, e in EVENT_IDS.items():
            label = f"Treat. $\\times$ 1(Year {e})"
            if label in st:
                for col, (b, se) in enumerate(zip(*st[label]), 1):
                    q = py[(py.table == "evi_median_mean") & (py.col == colmap[col]) & (py.term == f"event{e:+d}")]
                    if len(q): checks.append(("evi_median_mean", colmap[col], f"event{e:+d}", b, se, q.b.iloc[0], q.se.iloc[0]))
    if not checks:
        print("\n(no Stata tables found in output/tables - run code/run_demo.do first to compare)")
        return
    cmp = pd.DataFrame(checks, columns=["table", "col", "term", "stata_b", "stata_t", "python_b", "python_se"])
    cmp["python_t"] = cmp.python_b / cmp.python_se
    cmp["diff_b"] = (cmp.python_b - cmp.stata_b).round(3)
    cmp["diff_t"] = (cmp.python_t - cmp.stata_t).round(2)
    cmp = cmp[["table", "col", "term", "stata_b", "python_b", "diff_b", "stata_t", "python_t", "diff_t"]]
    print("\n" + "=" * 100 + "\nCOMPARISON WITH STATA (coefficients and t statistics read from output/tables/*.tex;\n"
          "esttab prints b to 3 decimals and t to 2 decimals, so differences within rounding are expected)\n" + "=" * 100)
    print(cmp.round(3).to_string(index=False))
    exact = cmp[~cmp.col.astype(str).str.startswith("ETWFE")]
    print(f"\nexactly comparable specs (OLS, PPML, Sun & Abraham): max |diff b| = {exact.diff_b.abs().max():.3f}, "
          f"max |diff t| = {exact.diff_t.abs().max():.2f}")
    print("ETWFE columns: jwdid interacts the GCP covariate with cohort/period dummies; the Python ETWFE enters it "
          "linearly, so small differences are expected there. Sun & Abraham t statistics: eventstudyinteract adds the "
          "sampling variance of the cohort shares, which the Python version omits.")
    cmp.to_csv(OUT / "python_vs_stata.csv", index=False)


if __name__ == "__main__":
    main()
