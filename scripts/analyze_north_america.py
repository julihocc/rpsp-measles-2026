"""
Annual measles cases, population and MCV1 coverage for Canada, Mexico and the United States,
1980-2024, from the archived Our World in Data (OWID) files in data/.

Usage:
    uv run python scripts/analyze_north_america.py

Every input is checked against its checksum manifest in data/ before anything is computed.
Outputs (under results/):
- figures/fig_cases_per100k.pdf: annual reported cases per 100,000 population, 1980-2024.
- tables/association.tex: Spearman correlation between annual cases and MCV1 coverage, with
  moving-block bootstrap intervals.
- tables/robustness.tex: the same correlation for year-on-year changes and longer blocks.
- results.json: rates per 100,000 by period, the correlations, the smallest correlation the
  tests could detect, and the contrast between Mexico's periods before and after 1992.
"""

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
COUNTRIES = ["Canada", "Mexico", "United States"]
FIRST_YEAR, LAST_YEAR = 1980, 2024
SEED = 20261101
PERIODS = [(1980, 1989), (1990, 1999), (2000, 2009), (2010, 2019), (2020, 2024)]

sns.set_style("whitegrid")


def snapshot_file(prefix: str, snapshot: str, data_dir: Path = DATA) -> Path:
    if not re.fullmatch(r"\d{8}", snapshot):
        raise ValueError("Snapshot must be an eight-digit date in YYYYMMDD format")
    path = data_dir / f"{prefix}_{snapshot}.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing input for snapshot {snapshot}: {path}")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_snapshot(manifest_path: Path, paths: list[Path]) -> None:
    """Check names, sizes and SHA-256 hashes of the inputs against a preserved manifest."""
    if not manifest_path.is_file():
        raise FileNotFoundError(f"No preserved checksum manifest: {manifest_path}")
    manifest = pd.read_csv(manifest_path)
    expected = manifest.set_index("file")[["bytes", "sha256"]]
    if set(expected.index) != {path.name for path in paths}:
        raise ValueError(f"{manifest_path} does not list exactly the inputs {sorted(p.name for p in paths)}")
    for path in paths:
        row = expected.loc[path.name]
        if path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            raise ValueError(f"Input does not match preserved checksum manifest: {path}")


def load(path: Path, value_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["Entity"].isin(COUNTRIES)].rename(columns={value_col: "value"})
    return df[["Entity", "Year", "value"]]


def block_starts(years, block: int = 3) -> np.ndarray:
    """Positions where a block of `block` consecutive calendar years starts."""
    years = np.asarray(years)
    return np.array(
        [i for i in range(len(years) - block + 1) if np.all(np.diff(years[i : i + block]) == 1)]
    )


def block_bootstrap_spearman(years, x, y, block=3, n_boot=5000, rng=None):
    """Bootstrap from length-three blocks containing consecutive calendar years."""
    years = np.asarray(years)
    n = len(x)
    starts = block_starts(years, block)
    if not len(starts):
        raise ValueError(f"No contiguous calendar-year blocks of length {block}")
    if rng is None:
        rng = np.random.default_rng()
    k = int(np.ceil(n / block))
    out = []
    for _ in range(n_boot):
        idx = np.concatenate([np.arange(s, s + block) for s in rng.choice(starts, k)])[:n]
        r = stats.spearmanr(x[idx], y[idx])[0]
        if not np.isnan(r):
            out.append(r)
    return np.percentile(out, [2.5, 97.5])


def tex_num(v: float) -> str:
    """Format a number for LaTeX so negative values get a true minus sign (and -0.00 prints as 0.00)."""
    return f"${round(v, 2) + 0.0:.2f}$"


def tex_ci(lo: float, hi: float) -> str:
    return f"$[{round(lo, 2) + 0.0:.2f}, {round(hi, 2) + 0.0:.2f}]$"


def extremes(s: pd.Series, how: str) -> str:
    """Extreme value of a series with every year it occurs (ties listed)."""
    v = s.min() if how == "min" else s.max()
    years = [int(y) for y in s.index[s == v]]
    if len(years) > 3:
        return f"{v:.0f} ({len(years)} years, latest {max(years)})"
    return f"{v:.0f} ({', '.join(map(str, years))})"


def missing_years(years: set, span: set) -> str:
    """Years absent from a series' span, summarised when there are many."""
    gaps = sorted(span - years)
    if not gaps:
        return "none"
    if len(gaps) > 5:
        return f"{len(gaps)} years in {gaps[0]}--{gaps[-1]}"
    return ", ".join(map(str, gaps))


def spearman_row(years, x, y, rng, block=3) -> dict:
    r = float(stats.spearmanr(x, y)[0])
    lo, hi = block_bootstrap_spearman(years, x, y, block=block, rng=rng)
    years = np.asarray(years)
    starts = block_starts(years, block)
    covered = {int(years[i]) for s in starts for i in range(s, s + block)}
    return {"n": len(x), "rho": r, "ci_low": float(lo), "ci_high": float(hi),
            "n_block_starts": len(starts),
            "years_never_resampled": sorted(int(y) for y in years if int(y) not in covered)}


def min_detectable_spearman(n: int, alpha: float = 0.05, power: float = 0.80) -> float:
    """Smallest |rho| a two-sided Spearman test detects with the given power.

    Uses the Fisher z approximation with the Bonett-Wright (2000) variance for
    Spearman's rho, (1 + rho^2 / 2) / (n - 3), and assumes independent
    observations, so autocorrelated annual series have less power than this.
    """
    z = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    return float(optimize.brentq(
        lambda r: np.arctanh(r) - z * np.sqrt((1 + r**2 / 2) / (n - 3)), 1e-9, 1 - 1e-9))


def tex_table(df: pd.DataFrame, path: Path, caption: str, label: str, cols_fmt: str):
    body = df.to_latex(index=False, escape=False, column_format=cols_fmt, na_rep="--")
    path.write_text(
        "\\begin{table}[htbp]\n\\centering\n\\small\n"
        f"\\caption{{{caption}}}\n\\label{{{label}}}\n{body}\\end{{table}}\n",
        encoding="utf-8",
        newline="\n",
    )


def load_inputs(data_dir: Path, snapshot: str, pop_snapshot: str) -> dict:
    """Validate the case, coverage and population files against their manifests and load them."""
    who_p = snapshot_file("reported_cases_measles", snapshot, data_dir)
    mcv1_p = snapshot_file("measles_vaccination_mcv1", snapshot, data_dir)
    mcv2_p = snapshot_file("measles_vaccination_mcv2", snapshot, data_dir)
    usl_p = snapshot_file("measles_us_cases", snapshot, data_dir)
    validate_snapshot(data_dir / f"data_provenance_{snapshot}.csv", [who_p, mcv1_p, mcv2_p, usl_p])
    pop_p = snapshot_file("population", pop_snapshot, data_dir)
    proj_p = snapshot_file("population_un_projections", pop_snapshot, data_dir)
    validate_snapshot(data_dir / f"data_provenance_population_{pop_snapshot}.csv", [pop_p, proj_p])
    usl = pd.read_csv(usl_p)
    usl = usl[["Year", "Number of measles cases"]].rename(columns={"Number of measles cases": "cdc"})
    annotation = pd.read_csv(usl_p)["Number of measles cases (Annotations)"].dropna().unique().tolist()
    return {
        "who": load(who_p, "Measles - number of reported cases"),
        "mcv1": load(mcv1_p, "Measles, first dose (MCV1)"),
        "mcv2": load(mcv2_p, "Measles, second dose (MCV2)"),
        "usl": usl,
        "annotation": annotation,
        "paths": [who_p, mcv1_p, mcv2_p, usl_p, pop_p, proj_p],
        "pop_p": pop_p,
        "proj_p": proj_p,
        "snapshot": snapshot,
        "pop_snapshot": pop_snapshot,
    }


def frames(inp: dict) -> dict:
    """Cases, population, rates and coverage on the 1980-2024 year index."""
    cases = inp["who"].pivot(index="Year", columns="Entity", values="value")
    cases = cases.loc[FIRST_YEAR:LAST_YEAR]
    # The US (CDC-sourced) series on the same year index, for sensitivity analyses
    cdc = inp["usl"].set_index("Year")["cdc"].reindex(cases.index)

    # Population from OWID (UN WPP 2024 estimates up to 2023; 2024 is the WPP 2024 medium-variant
    # projection, the only year the estimates do not cover)
    pop = load(inp["pop_p"], "Population").pivot(index="Year", columns="Entity", values="value")
    proj = pd.read_csv(inp["proj_p"])
    proj = proj[proj.Entity.isin(COUNTRIES)]
    est = proj.pivot(index="Year", columns="Entity", values="Population").loc[FIRST_YEAR:LAST_YEAR]
    common = est.dropna().index
    if not np.array_equal(est.loc[common, COUNTRIES].values, pop.loc[common, COUNTRIES].values):
        raise ValueError("The two population files disagree for the years they share")
    projected = proj.pivot(index="Year", columns="Entity", values="Population (Projected)")
    popw = pop.reindex(cases.index)
    projected_years = [int(y) for y in popw.index[popw[COUNTRIES].isna().any(axis=1)]]
    popw = popw.fillna(projected.reindex(cases.index))
    if popw[COUNTRIES].isna().any().any():
        raise ValueError("Population missing for some years of the analysis window")
    return {
        "cases": cases,
        "cdc": cdc,
        "popw": popw,
        "projected_years": projected_years,
        "rate": cases[COUNTRIES] / popw[COUNTRIES] * 1e5,
        "rate_cdc": cdc / popw["United States"] * 1e5,
        "v1": inp["mcv1"].pivot(index="Year", columns="Entity", values="value"),
        "v2": inp["mcv2"].pivot(index="Year", columns="Entity", values="value"),
        "case_series": [(c, cases[c]) for c in COUNTRIES] + [("United States (US series)", cdc)],
    }


def shared_results(inp: dict, fr: dict) -> dict:
    """Rates by period, the coverage-cases correlations, their power and Mexico's two periods."""
    access_date = re.search(r"_(\d{8})\.csv", inp["paths"][0].name).group(1)
    res = {"access_date": access_date, "population_access_date": inp["pop_snapshot"], "seed": SEED,
           "us_series_annotation": inp["annotation"]}

    rate, rate_cdc = fr["rate"], fr["rate_cdc"]
    pc_rows = []
    for c, full in [(c, rate[c]) for c in COUNTRIES] + [("United States (US series)", rate_cdc)]:
        s = full.dropna()
        for d0, d1 in PERIODS:
            sub = s[(s.index >= d0) & (s.index <= d1)]
            pc_rows.append({"Country": c, "Period": f"{d0}--{d1}", "n years": len(sub),
                            "Median": float(sub.median()), "Max": float(sub.max()),
                            "Year of max": int(sub.idxmax())})
    res["per_capita"] = {
        "population_source": "OWID population (UN WPP 2024 estimates); WPP 2024 medium-variant projection "
                             "for years without an estimate",
        "projected_years": fr["projected_years"],
        "decades_per100k": pc_rows,
    }

    # --- Association between MCV1 and cases (ecological, descriptive)
    # Spearman is rank-based, so raw counts and log counts give identical results.
    # Rows are drawn in a fixed order from one generator; new rows go last so earlier rows keep their intervals.
    v1, case_series = fr["v1"], fr["case_series"]

    def window(d: pd.DataFrame) -> pd.DataFrame:
        return d[(d.index >= FIRST_YEAR) & (d.index <= LAST_YEAR)]

    rng = np.random.default_rng(SEED)
    arows = []
    for c, y in case_series:
        mcv = v1["United States" if c.startswith("United States") else c]
        for lag in (0, 1):
            d = window(pd.DataFrame({"y": y, "x": mcv.shift(lag)}).dropna())
            arows.append({"Country": c, "MCV1 lag (y)": lag,
                          **spearman_row(d.index, d.x.values, d.y.values, rng)})
    res["association"] = arows

    # --- Robustness: shared long-run trends and block length (separate generator)
    rng_rob = np.random.default_rng(SEED + 1)
    rrows = []
    for c, y in case_series:
        mcv = v1["United States" if c.startswith("United States") else c]
        # Year-on-year changes; differences next to a missing year are dropped
        d = window(pd.DataFrame({"y": np.log10(1 + y).diff(), "x": mcv.diff()}).dropna())
        rrows.append({"Country": c, "Analysis": "First differences, lag 0",
                      **spearman_row(d.index, d.x.values, d.y.values, rng_rob)})
        d = window(pd.DataFrame({"y": y, "x": mcv}).dropna())
        for block in (5, 8):
            rrows.append({"Country": c, "Analysis": f"Levels, lag 0, block length {block}",
                          **spearman_row(d.index, d.x.values, d.y.values, rng_rob, block=block)})
    # Added later, so drawn after the rows above (which keep their intervals): change in cases
    # against the previous year's change in MCV1, the first-difference version of the lag-1 rows
    # in the association table.
    for c, y in case_series:
        mcv = v1["United States" if c.startswith("United States") else c]
        d = window(pd.DataFrame({"y": np.log10(1 + y).diff(), "x": mcv.shift(1).diff()}).dropna())
        rrows.append({"Country": c, "Analysis": "First differences, lag 1",
                      **spearman_row(d.index, d.x.values, d.y.values, rng_rob)})
    # Display grouped by country, first differences before levels (stable sort keeps lag/block order)
    country_order = [c for c, _ in case_series]
    rrows.sort(key=lambda r: (country_order.index(r["Country"]), r["Analysis"].startswith("Levels")))
    res["robustness"] = rrows

    # Power: smallest correlation the tests above could detect (closed form, draws no random numbers)
    res["power"] = {
        "method": "Fisher z with Bonett-Wright variance for Spearman rho, independent observations",
        "alpha_two_sided": 0.05, "power": 0.80,
        "min_detectable_abs_rho_by_n": {int(n): round(min_detectable_spearman(int(n)), 4)
                                        for n in sorted({r["n"] for r in arows + rrows})},
        "illustration": {"n": 40, "min_detectable_abs_rho": round(min_detectable_spearman(40), 4)},
    }

    # Mexico's levels correlation contrasts two periods: low coverage and many cases up to 1991,
    # high coverage and few cases from 1992 (the first year with fewer than 1000 cases).
    mx = window(pd.DataFrame({"cases": fr["cases"]["Mexico"], "mcv1": v1["Mexico"]}).dropna())
    low = mx.index[mx.cases < 1000]
    split = int(low.min()) if len(low) and low.min() > mx.index.min() else None
    res["mexico_periods"] = {} if split is None else {
        f"{lo}-{hi}": {"mcv1_min": float(p.mcv1.min()), "mcv1_max": float(p.mcv1.max()),
                       "cases_min": int(p.cases.min()), "cases_max": int(p.cases.max()),
                       "cases_median": float(p.cases.median()), "n_years": len(p)}
        for lo, hi in [(int(mx.index.min()), split - 1), (split, int(mx.index.max()))]
        for p in [mx.loc[lo:hi]]
    }
    return res


def write_shared_outputs(res: dict, fr: dict, out_dir: Path) -> None:
    """Write the association and robustness tables and the per-capita figure."""
    tab, fig_dir = out_dir / "tables", out_dir / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    tex_table(pd.DataFrame([{"Country": r["Country"], "MCV1 lag (y)": r["MCV1 lag (y)"], "n": r["n"],
                             "Spearman rho": tex_num(r["rho"]),
                             "95\\% block-bootstrap CI": tex_ci(r["ci_low"], r["ci_high"])}
                            for r in res["association"]]),
              tab / "association.tex",
              "Spearman correlation between annual reported cases and MCV1 coverage, "
              "with a 95\\% moving-block bootstrap interval (blocks of three consecutive years, 5000 resamples). "
              "WHO series for each country; the last two rows use the OWID US series. Ecological and descriptive.",
              "tab:assoc", "lrrrr")
    tex_table(pd.DataFrame([{"Country": r["Country"], "Analysis": r["Analysis"], "n": r["n"],
                             "Spearman rho": tex_num(r["rho"]),
                             "95\\% block-bootstrap CI": tex_ci(r["ci_low"], r["ci_high"])}
                            for r in res["robustness"]]),
              tab / "robustness.tex",
              "Robustness of the association in Table~\\ref{tab:assoc}. First differences, which remove shared "
              "long-run trends: year-on-year change in $\\log_{10}(1+\\text{cases})$ versus the change in MCV1 "
              "in the same year (lag 0) or the previous year (lag 1), block length 3. Levels: as in "
              "Table~\\ref{tab:assoc} with longer bootstrap blocks.",
              "tab:robust", "llrrr")

    # One is added to cases before dividing, so years with zero cases stay on the log scale
    cases, cdc, popw = fr["cases"], fr["cdc"], fr["popw"]
    fig, ax = plt.subplots(figsize=(12, 6))
    for c in COUNTRIES:
        s = cases[c].dropna()
        ax.plot(s.index, (s.values + 1) / popw.loc[s.index, c].values * 1e5, marker="o", ms=3,
                label=f"{c} (WHO series)")
    s = cdc.dropna()
    ax.plot(s.index, (s.values + 1) / popw.loc[s.index, "United States"].values * 1e5, ls="--", color="k",
            lw=1.2, label="United States (OWID US series)")
    ax.set_yscale("log")
    ax.set_xlabel("Year")
    ax.set_ylabel("(Reported cases + 1) per 100 000 (log scale)")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, fig_dir / "fig_cases_per100k.pdf")


def save_figure(fig, path: Path, **kwargs) -> None:
    """Save a PDF without a creation date, so reruns give identical files, and close it."""
    fig.savefig(path, metadata={"CreationDate": None}, **kwargs)
    plt.close(fig)


def main(snapshot: str = "20261101", pop_snapshot: str = "20261107", data_dir: Path = DATA,
         out_dir: Path = RESULTS) -> dict:
    inp = load_inputs(data_dir, snapshot, pop_snapshot)
    fr = frames(inp)
    res = shared_results(inp, fr)
    write_shared_outputs(res, fr, out_dir)
    (out_dir / "results.json").write_text(json.dumps(res, indent=2, default=str), encoding="utf-8",
                                          newline="\n")
    print(pd.DataFrame(res["association"]).to_string())
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--snapshot", default="20261101", help="Snapshot date of the OWID case and coverage files")
    parser.add_argument("--population-snapshot", default="20261107",
                        help="Snapshot date of the two OWID population files")
    parser.add_argument("--data-dir", type=Path, default=DATA, help="Folder with the input files and manifests")
    parser.add_argument("--out-dir", type=Path, default=RESULTS, help="Folder for the outputs")
    args = parser.parse_args()
    main(args.snapshot, args.population_snapshot, args.data_dir, args.out_dir)
