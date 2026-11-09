"""
Monthly analysis of the 2025-2026 measles resurgence in Canada, Mexico and the United States.

Inputs are the archived, checksummed files in data/; nothing is downloaded:
- WHO provisional monthly measles file (who_measles_monthly_20261107.xlsx): monthly cases.
- WHO annual series (who_measles_annual_20261107.csv): 2025 totals for Canada and Mexico.
- OWID reported cases (WHO series) and US series (CDC), snapshot 20261101: 1980s and 2010s
  rates, and the US 2025 total.
- OWID population and population with UN projections, snapshot 20261107.
- OWID MCV1 and MCV2 coverage, snapshot 20261101.

Method, fixed before looking at the fits:
- Analysis window October 2024 - September 2026. October 2024 is when Canada's outbreak began and
  the first month after Canada's last missing monthly report (September 2024); nothing is imputed.
- Wave phases follow a mechanical rule. A peak is a month with at least 100 cases that is the
  maximum of the months within three months either side, and is not the last month of the window.
  The trough before the first peak is the lowest month from the window start to that peak; the
  trough between two peaks is the lowest month between them; after the last peak it is the lowest
  month to the window end. Ties go to the later month. Growth runs from a trough to the next peak,
  decline from a peak to the next trough. If the window ends above the last trough, the months
  from that trough to the end form an ongoing growth phase.
- Each phase is fitted with a quasi-Poisson log-linear model of monthly cases on the mid-month
  date, with log(days in month) as offset. The slope r is the daily growth rate; 95% intervals use
  the Pearson dispersion and a t distribution with n - 2 degrees of freedom.
- r is converted to a reproduction number per generation, R = exp(r * Tg), for a fixed generation
  interval Tg of 10, 12 and 14 days (Wallinga and Lipsitch 2007; mean measles serial interval about
  12 days, Vink et al. 2014).
- The file is provisional (reported to WHO as of October 2026), so every phase is refitted with the
  window ending in July 2026, leaving out August and September 2026.

Outputs under results/: tables/resurgence_*.tex, figures/fig_monthly.pdf and
results_resurgence.json. Run scripts/analyze_north_america.py first: the 1980s and 2010s rates
are checked against its results/results.json.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_north_america as ana  # noqa: E402
import check_who_series as cws  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
ISO = {"CAN": "Canada", "MEX": "Mexico", "USA": "United States"}
WINDOW = ("2024-10-31", "2026-09-30")
SENSITIVITY_END = "2026-07-31"
PEAK_MIN, HALF_WIDTH = 100, 3
TG_DAYS = (10, 12, 14)
CASES_SNAPSHOT, POP_SNAPSHOT, WHO_SNAPSHOT = "20261101", "20261107", "20261107"


# --- Data -----------------------------------------------------------------------------------

def monthly_series(path: Path) -> pd.DataFrame:
    """Monthly 'Measles total' (clinical, epi-linked and lab-confirmed) by country, month-end index."""
    m = pd.read_excel(path, sheet_name="WEB")
    m.columns = [" ".join(str(c).split()) for c in m.columns]
    m = m[m["ISO3"].isin(ISO)]
    m = m.assign(date=pd.to_datetime(dict(year=m.Year, month=m.Month, day=1)) + pd.offsets.MonthEnd(0))
    wide = m.pivot_table(index="date", columns="ISO3", values="Measles total", aggfunc="sum")
    full = pd.date_range(wide.index.min(), wide.index.max(), freq="ME")
    return wide.reindex(full)[list(ISO)]


def population(pop_p: Path, proj_p: Path) -> pd.DataFrame:
    """Population by year: UN WPP 2024 estimates, then the medium-variant projection."""
    est = pd.read_csv(pop_p)
    est = est[est.Entity.isin(ISO.values())].pivot(index="Year", columns="Entity", values="Population")
    proj = pd.read_csv(proj_p)
    proj = proj[proj.Entity.isin(ISO.values())]
    p_est = proj.pivot(index="Year", columns="Entity", values="Population")
    p_prj = proj.pivot(index="Year", columns="Entity", values="Population (Projected)")
    common = est.dropna().index.intersection(p_est.dropna().index)
    if not np.array_equal(est.loc[common].values, p_est.loc[common, est.columns].values):
        raise ValueError("The two population files disagree for the years they share")
    return est.combine_first(p_prj)


# --- Wave phases ----------------------------------------------------------------------------

def last_argmin(v: np.ndarray, a: int, b: int) -> int:
    """Index of the lowest value in v[a..b] (inclusive); ties go to the later month."""
    seg = v[a : b + 1]
    return a + int(np.flatnonzero(seg == seg.min())[-1])


def find_peaks(v: np.ndarray) -> list[int]:
    peaks = []
    for i in range(len(v) - 1):  # the last month cannot be a peak: the next month is unknown
        lo, hi = max(0, i - HALF_WIDTH), min(len(v), i + HALF_WIDTH + 1)
        if v[i] >= PEAK_MIN and v[i] == v[lo:hi].max() and not (peaks and i - peaks[-1] <= HALF_WIDTH):
            peaks.append(i)
    return peaks


def find_phases(v: np.ndarray) -> list[dict]:
    """Growth and decline phases (start and end positions, inclusive) from the rule in the docstring."""
    n, peaks = len(v), find_peaks(v)
    if not peaks:
        return []
    troughs = [last_argmin(v, 0, peaks[0])]
    troughs += [last_argmin(v, p + 1, q - 1) for p, q in zip(peaks, peaks[1:])]
    troughs.append(last_argmin(v, peaks[-1] + 1, n - 1))
    phases = []
    for k, p in enumerate(peaks):
        phases.append({"wave": k + 1, "phase": "growth", "i0": troughs[k], "i1": p})
        phases.append({"wave": k + 1, "phase": "decline", "i0": p, "i1": troughs[k + 1]})
    if troughs[-1] < n - 1 and v[-1] > v[troughs[-1]]:
        phases.append({"wave": len(peaks) + 1, "phase": "growth (ongoing)", "i0": troughs[-1], "i1": n - 1})
    return phases


def mid_month(index: pd.DatetimeIndex) -> pd.DatetimeIndex:
    """Middle of each month of a month-end index."""
    return index - pd.to_timedelta(index.days_in_month.to_numpy(dtype=float) / 2, unit="D")


def glm_fit(s: pd.Series):
    """Quasi-Poisson GLM of monthly cases on days since the first mid-month, offset log(days in month)."""
    days = s.index.days_in_month.to_numpy(dtype=float)
    mid = mid_month(s.index).to_numpy()
    t = (mid - mid[0]) / np.timedelta64(1, "D")
    model = sm.GLM(s.to_numpy(dtype=float), sm.add_constant(t), family=sm.families.Poisson(),
                   offset=np.log(days))
    return model.fit(scale="X2")


def fit_phase(s: pd.Series) -> dict:
    """Daily growth rate of one phase, its 95% interval and the implied reproduction numbers."""
    y = s.to_numpy(dtype=float)
    n = len(y)
    out = {"n_months": n, "first": str(s.index[0].date())[:7], "last": str(s.index[-1].date())[:7],
           "cases_first": int(y[0]), "cases_last": int(y[-1])}
    if n < 3:
        out["month_ratio"] = float((y[-1] + 1) / (y[0] + 1))
        return out
    res = glm_fit(s)
    r, se = float(res.params[1]), float(res.bse[1])
    q = stats.t.ppf(0.975, n - 2)
    out.update({"r_per_day": r, "r_ci": [r - q * se, r + q * se], "dispersion": float(res.scale),
                "R_by_Tg": {str(tg): float(np.exp(r * tg)) for tg in TG_DAYS},
                "R_Tg12_ci": [float(np.exp((r - q * se) * 12)), float(np.exp((r + q * se) * 12))]})
    return out


def phases_for(series: dict[str, pd.Series], end: str) -> list[dict]:
    rows = []
    for iso, s in series.items():
        s = s.loc[:end]
        for ph in find_phases(s.to_numpy(dtype=float)):
            fit = fit_phase(s.iloc[ph["i0"] : ph["i1"] + 1])
            rows.append({"country": ISO[iso], "iso": iso, "wave": ph["wave"], "phase": ph["phase"],
                         "i0": ph["i0"], "i1": ph["i1"], **fit})
    return rows


# --- Output helpers -------------------------------------------------------------------------

def num(v: float) -> str:
    return f"{int(round(v)):,}"


def signed(v: float) -> str:
    """Rounded signed number for LaTeX text, with a typeset minus sign."""
    s = f"{v:+.0f}"
    return "$-$" + s[1:] if s.startswith("-") else s


def pct_week(r: float) -> float:
    return 100 * (np.exp(7 * r) - 1)


def write_table(path: Path, caption: str, label: str, header: list[str], rows: list[list[str]],
                cols: str, note: str | None = None, wide: bool = False) -> None:
    """LaTeX table; `wide` spans the full page width (adjustwidth with \\extralength)."""
    lines = ["\\begin{table}[H]"]
    if wide:
        lines += ["\\begin{adjustwidth}{-\\extralength}{0cm}", "\\setlength{\\tabcolsep}{4pt}"]
    lines += [f"\\caption{{{caption}\\label{{{label}}}}}", "\\centering", "\\small",
              f"\\begin{{tabular}}{{{cols}}}", "\\toprule",
              " & ".join(f"\\textbf{{{h}}}" for h in header) + "\\\\", "\\midrule"]
    lines += [" & ".join(r) + "\\\\" for r in rows]
    lines += ["\\bottomrule", "\\end{tabular}"]
    if note:
        lines += ["\\par\\smallskip", "{\\raggedright\\footnotesize " + note + "\\par}"]
    if wide:
        lines += ["\\end{adjustwidth}"]
    lines += ["\\end{table}", ""]
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def phase_label(month: str) -> str:
    return pd.Timestamp(month + "-01").strftime("%b %Y")


# --- Main -----------------------------------------------------------------------------------

def main(data_dir: Path = DATA, out_dir: Path = RESULTS) -> None:
    descriptive = out_dir / "results.json"
    if not descriptive.is_file():
        raise SystemExit(f"Run scripts/analyze_north_america.py first: {descriptive} not found")
    FIG, TAB = out_dir / "figures", out_dir / "tables"
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    monthly_p = data_dir / f"who_measles_monthly_{WHO_SNAPSHOT}.xlsx"
    annual_p = data_dir / f"who_measles_annual_{WHO_SNAPSHOT}.csv"
    cws.validate(data_dir / f"data_provenance_who_{WHO_SNAPSHOT}.csv", [monthly_p, annual_p])
    case_paths = {k: ana.snapshot_file(k, CASES_SNAPSHOT, data_dir) for k in
                  ["reported_cases_measles", "measles_us_cases", "measles_vaccination_mcv1",
                   "measles_vaccination_mcv2"]}
    cws.validate(data_dir / f"data_provenance_{CASES_SNAPSHOT}.csv", list(case_paths.values()), exact=False)
    pop_p = ana.snapshot_file("population", POP_SNAPSHOT, data_dir)
    proj_p = ana.snapshot_file("population_un_projections", POP_SNAPSHOT, data_dir)
    cws.validate(data_dir / f"data_provenance_population_{POP_SNAPSHOT}.csv", [pop_p, proj_p])

    monthly = monthly_series(monthly_p)
    as_of = cws.monthly_as_of(monthly_p)
    pop = population(pop_p, proj_p)
    who = ana.load(case_paths["reported_cases_measles"], "Measles - number of reported cases")
    who = who.pivot(index="Year", columns="Entity", values="value")
    usl = pd.read_csv(case_paths["measles_us_cases"]).set_index("Year")["Number of measles cases"]
    annual = cws.load_annual(annual_p).set_index(["iso3", "year"])["who_annual"]

    res = {"monthly_file_as_of": as_of, "window": list(WINDOW), "sensitivity_end": SENSITIVITY_END,
           "peak_rule": {"min_cases": PEAK_MIN, "half_width_months": HALF_WIDTH}, "Tg_days": list(TG_DAYS)}

    # Monthly series in the window; nothing is imputed
    win = monthly.loc[WINDOW[0] : WINDOW[1]]
    missing = {ISO[c]: [str(d.date())[:7] for d in win.index[win[c].isna()]] for c in ISO}
    if any(missing.values()):
        raise ValueError(f"Months missing inside the analysis window: {missing}")
    series = {c: win[c] for c in ISO}
    res["monthly_cases"] = {ISO[c]: {str(d.date())[:7]: int(v) for d, v in s.items()} for c, s in series.items()}
    res["missing_months_before_window"] = {
        ISO[c]: [str(d.date())[:7] for d in monthly.loc[:WINDOW[0]].index[monthly.loc[:WINDOW[0], c].isna()]]
        for c in ISO}

    # Reported rates per 100,000: 1980s and 2010s medians against 2025 and January-September 2026
    jan_sep = monthly.loc["2026-01-31":"2026-09-30"]
    if jan_sep.isna().any().any() or len(jan_sep) != 9:
        raise ValueError("January-September 2026 is incomplete in the monthly file")
    rate_rows, tex_rows = [], []
    for iso, name in ISO.items():
        hist = usl if iso == "USA" else who[name]  # US: CDC-sourced series throughout
        r = hist.loc[1980:2019] / pop[name].loc[1980:2019] * 1e5
        c2025 = float(usl.loc[2025]) if iso == "USA" else float(annual.loc[(iso, 2025)])
        c2026 = float(jan_sep[iso].sum())
        row = {"country": name, "source_history_and_2025": "US series (CDC)" if iso == "USA" else "WHO series",
               "median_1980s": float(r.loc[1980:1989].median()), "median_2010s": float(r.loc[2010:2019].median()),
               "cases_2025": int(c2025), "rate_2025": c2025 / pop[name].loc[2025] * 1e5,
               "cases_2026_jan_sep": int(c2026), "rate_2026_jan_sep": c2026 / pop[name].loc[2026] * 1e5}
        if iso == "USA":
            row["cdc_2026_to_date"] = int(usl.loc[2026])
            row["rate_cdc_2026_to_date"] = float(usl.loc[2026]) / pop[name].loc[2026] * 1e5
        # Last earlier year with an annual rate at least as high (a nine-month 2026 total is compared
        # with full years, so "higher than every year since" is conservative for 2026)
        full = (hist.loc[:2024] / pop[name].reindex(hist.loc[:2024].index) * 1e5).dropna()
        for key in ["rate_2025", "rate_2026_jan_sep"] + (["rate_cdc_2026_to_date"] if iso == "USA" else []):
            above = full[full >= row[key]]
            row[f"last_year_at_or_above_{key}"] = int(above.index.max()) if len(above) else None
        row["history_first_year"] = int(full.index.min())
        rate_rows.append(row)
        tex_rows.append([name, f"{row['median_1980s']:.3g}", f"{row['median_2010s']:.2g}",
                         f"{num(c2025)} ({row['rate_2025']:.3g})", f"{num(c2026)} ({row['rate_2026_jan_sep']:.3g})"])
    # Check against the descriptive analysis (same files, same definitions)
    ref = {(d["Country"], d["Period"]): d["Median"] for d in
           json.loads(descriptive.read_text())["per_capita"]["decades_per100k"]}
    for row in rate_rows:
        key = "United States (US series)" if row["country"] == "United States" else row["country"]
        for per, col in [("1980--1989", "median_1980s"), ("2010--2019", "median_2010s")]:
            if not np.isclose(ref[(key, per)], row[col], rtol=1e-9, atol=1e-12):
                raise ValueError(f"Rate for {key} {per} differs from {descriptive}")
    res["rates_per100k"] = rate_rows
    us = rate_rows[-1]
    write_table(
        TAB / "resurgence_rates.tex",
        "Reported measles cases per 100,000 population: median annual rate in 1980--1989 and 2010--2019, and "
        "cases (rate) in 2025 and in January--September 2026.",
        "tab:rates",
        ["Country", "1980--1989", "2010--2019", "2025", "Jan--Sep 2026"],
        tex_rows, "lrrrr",
        note=("Cases: WHO series (Canada, Mexico) and CDC-sourced US series (United States) to 2025; WHO "
              f"provisional monthly data for 2026 (as of {as_of}). The 2026 column covers nine months. The "
              f"CDC total for 2026 to date is {num(us['cdc_2026_to_date'])} ({us['rate_cdc_2026_to_date']:.3g} per "
              "100,000). Population: UN World Population Prospects 2024, medium-variant projection for 2025 and 2026."),
    )

    # Wave phases, full window and window ending July 2026
    main_rows = phases_for(series, WINDOW[1])
    sens_rows = phases_for(series, SENSITIVITY_END)
    res["phases"], res["phases_to_july_2026"] = main_rows, sens_rows
    sens = {(r["country"], r["wave"], r["phase"]): r for r in sens_rows}
    tex_rows = []
    for r in main_rows:
        s = sens.get((r["country"], r["wave"], r["phase"]))
        months = f"{phase_label(r['first'])}--{phase_label(r['last'])}"
        if "r_per_day" not in r:
            raise ValueError(f"Phase with fewer than three months: {r}")
        g = pct_week(r["r_per_day"])
        lo, hi = (pct_week(x) for x in r["r_ci"])
        t2 = np.log(2) / abs(r["r_per_day"]) / 7
        Rs = [r["R_by_Tg"][str(tg)] for tg in TG_DAYS]
        if s is None:
            s_txt = "--"
        elif s["first"] == r["first"] and s["last"] == r["last"]:
            s_txt = "same"
        else:
            s_txt = f"{signed(pct_week(s['r_per_day']))} ({phase_label(s['last'])})"
        first_of_country = r is main_rows[0] or r["country"] != prev
        tex_rows.append([r["country"] if first_of_country else "",
                         f"{r['wave']}, {r['phase']}", months, str(r["n_months"]),
                         f"{signed(g)} ({signed(lo)}, {signed(hi)})", f"{t2:.1f}", f"{min(Rs):.2f}--{max(Rs):.2f}",
                         s_txt])
        prev = r["country"]
    write_table(
        TAB / "resurgence_phases.tex",
        "Growth and decline phases of reported measles cases by month, October 2024--September 2026.",
        "tab:phases",
        ["Country", "Wave, phase", "Months", "n", "\\% per week (95\\% CI)", "$T_2$ (weeks)",
         "$R$ ($T_g$ 10--14 d)", "To Jul 2026"],
        tex_rows, "lllrrrrr", wide=True,
        note=("Phases follow the mechanical rule in Section~\\ref{sec:methods}. \\% per week: change in "
              "monthly cases per week, $100(e^{7r}-1)$, from a quasi-Poisson log-linear fit with the days in the "
              "month as offset; $T_2$: doubling (growth) or halving (decline) time. $R$: reproduction number per "
              f"generation, $e^{{rT_g}}$, for a fixed generation interval $T_g$ of {TG_DAYS[0]}--{TG_DAYS[-1]} days. "
              "To Jul 2026: the same fit with August and September 2026 left out (``same'' if the phase ends "
              "before August 2026; otherwise \\% per week and the new last month)."),
    )

    # MCV1 and MCV2 coverage, 2019-2025
    cov = {}
    rows = []
    for key, dose, col in [("measles_vaccination_mcv1", "MCV1", "Measles, first dose (MCV1)"),
                           ("measles_vaccination_mcv2", "MCV2", "Measles, second dose (MCV2)")]:
        d = ana.load(case_paths[key], col).pivot(index="Year", columns="Entity", values="value").loc[2019:2025]
        cov[dose] = {c: {int(y): (None if pd.isna(v) else int(v)) for y, v in d[c].items()} for c in ISO.values()}
    years = list(range(2019, 2026))
    for c in ISO.values():
        for dose in ["MCV1", "MCV2"]:
            rows.append([c if dose == "MCV1" else "", dose] +
                        ["--" if cov[dose][c].get(y) is None else str(cov[dose][c][y]) for y in years])
    res["coverage_2019_2025"] = cov
    write_table(TAB / "resurgence_coverage.tex",
                "Measles vaccine coverage (\\% of children, WUENIC via OWID), 2019--2025.",
                "tab:cov", ["Country", "Dose"] + [str(y) for y in years], rows, "ll" + "r" * len(years))

    # Figure: monthly cases with the fitted phases
    fig, axes = plt.subplots(3, 1, figsize=(10, 6.2), sharex=True)
    for ax, (iso, s) in zip(axes, series.items()):
        ax.plot(mid_month(s.index), s + 1, "o-", color="0.3", ms=4, lw=1, label="Reported cases + 1")
        for r in [x for x in main_rows if x["iso"] == iso]:
            seg = s.iloc[r["i0"] : r["i1"] + 1]  # the same months as the fit in the table
            color = "tab:red" if r["phase"].startswith("growth") else "tab:blue"
            ax.plot(mid_month(seg.index), glm_fit(seg).fittedvalues + 1, "-", color=color, lw=2.2, alpha=0.8)
        ax.axvspan(pd.Timestamp("2026-08-01"), pd.Timestamp("2026-10-01"), color="0.85", zorder=0)
        ax.set_yscale("log")
        ax.set_ylabel("Cases + 1 (log scale)")
        ax.set_title(ISO[iso], loc="left", fontsize=11)
    axes[0].plot([], [], "-", color="tab:red", lw=2.2, label="Fitted growth phase")
    axes[0].plot([], [], "-", color="tab:blue", lw=2.2, label="Fitted decline phase")
    axes[0].fill_between([], [], color="0.85", label="Aug–Sep 2026 (least complete)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=4, fontsize=9, frameon=False)
    axes[-1].set_xlabel("Month")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    ana.save_figure(fig, FIG / "fig_monthly.pdf")

    (out_dir / "results_resurgence.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"rates": rate_rows}, indent=1))
    for r in main_rows:
        print(r["country"], r["wave"], r["phase"], r["first"], r["last"], r["n_months"],
              f"{pct_week(r['r_per_day']):+.1f}%/wk", {k: round(v, 2) for k, v in r["R_by_Tg"].items()},
              f"disp={r['dispersion']:.1f}")
    print("--- to July 2026")
    for r in sens_rows:
        print(r["country"], r["wave"], r["phase"], r["first"], r["last"], f"{pct_week(r['r_per_day']):+.1f}%/wk")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", type=Path, default=DATA, help="Folder with the input files and manifests")
    parser.add_argument("--out-dir", type=Path, default=RESULTS, help="Folder for the outputs")
    args = parser.parse_args()
    main(args.data_dir, args.out_dir)
