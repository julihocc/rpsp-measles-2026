"""
Check WHO's annual measles series against WHO's monthly data and national totals.

For Canada, Mexico and the United States, 2015-2024, compare three counts per year:

- WHO annual: official reported cases (GHO WHS3_62, WHO/UNICEF Joint Reporting Form),
  the series Our World in Data republishes.
- WHO monthly: the sum of WHO's provisional monthly "Measles total" for that year.
- National: CDC for the United States (the CDC-sourced OWID US series in the 20261101
  snapshot), PHAC for Canada and DGE for Mexico (data/national_measles_totals.csv,
  one cited source per value).

The WHO annual value "matches" a count when it is within one case of it.

Usage:
    uv run python scripts/check_who_series.py

Every input is checked against its checksum manifest in data/ before anything is written.
Writes results/who_series_check.csv and the LaTeX table results/tables/who_check.tex, and
prints a per-country summary.

It also tabulates the preliminary 2025 and partial-2026 totals after the analysis window
(tables/post2024.tex, post2024_totals.csv): WHO annual, WHO monthly, CDC via the OWID US series,
and DGE values for Mexico from data/national_measles_totals_2025_2026.csv (one cited source and
cut-off per value; checksum in data_provenance_national_2025_2026.csv). No PHAC total for
2025-2026 was recorded, so Canada has no national value there.
"""

import argparse
import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
NATIONAL = "national_measles_totals.csv"
NATIONAL_RECENT = "national_measles_totals_2025_2026.csv"
RECENT_YEARS = (2025, 2026)
CDC_FILE, CDC_SNAPSHOT = "measles_us_cases_20261101.csv", "20261101"
COUNTRIES = {"CAN": "Canada", "MEX": "Mexico", "USA": "United States"}
FIRST_YEAR, LAST_YEAR = 2015, 2024


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(manifest_path: Path, paths: list[Path], exact: bool = True) -> None:
    """Check names, sizes and SHA-256 hashes against a preserved manifest."""
    expected = pd.read_csv(manifest_path).set_index("file")
    names = {p.name for p in paths}
    if exact and set(expected.index) != names:
        raise ValueError(f"{manifest_path} does not list exactly {sorted(names)}")
    for path in paths:
        row = expected.loc[path.name]
        if path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            raise ValueError(f"Input does not match {manifest_path}: {path}")


def load_monthly(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="WEB", dtype=str)
    df.columns = [" ".join(str(c).split()) for c in df.columns]  # "Measles \ntotal" -> "Measles total"
    df = df[df["ISO3"].isin(COUNTRIES)].copy()
    df["year"] = df["Year"].astype(int)
    df["total"] = pd.to_numeric(df["Measles total"]).fillna(0).astype(int)
    return (
        df.groupby(["ISO3", "year"])
        .agg(who_monthly=("total", "sum"), months_reported=("Month", "nunique"),
             first_month=("Month", "min"), last_month=("Month", "max"))
        .reset_index()
        .rename(columns={"ISO3": "iso3"})
    )


def load_annual(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[(df["spatial_type"] == "COUNTRY") & df["iso3"].isin(COUNTRIES)]
    return df[["iso3", "year", "cases"]].rename(columns={"cases": "who_annual"})


def load_national(data_dir: Path) -> pd.DataFrame:
    cdc_path = data_dir / CDC_FILE
    cdc = pd.read_csv(cdc_path).rename(columns={"Year": "year", "Number of measles cases": "national"})
    cdc = cdc.assign(iso3="USA", national_source="CDC (OWID number-of-measles-cases, 20261101 snapshot)")
    other = pd.read_csv(data_dir / NATIONAL)[["iso3", "year", "cases", "source"]].rename(
        columns={"cases": "national", "source": "national_source"}
    )
    return pd.concat([cdc[["iso3", "year", "national", "national_source"]], other], ignore_index=True)


def within_one(a, b) -> bool:
    return pd.notna(a) and pd.notna(b) and abs(a - b) <= 1


def write_tex_table(df: pd.DataFrame, path: Path, monthly_as_of: str) -> None:
    """One row per year; WHO annual, WHO monthly and national totals for each country."""
    fmt = lambda v: "--" if pd.isna(v) else f"{int(v)}"
    national_name = {"CAN": "PHAC", "MEX": "DGE", "USA": "CDC"}
    order = ["USA", "CAN", "MEX"]
    lines = []
    for year, g in df.groupby("year"):
        g = g.set_index("iso3")
        cells = [fmt(g.loc[iso, c]) for iso in order for c in ("who_annual", "who_monthly", "national")]
        lines.append(f"{year} & " + " & ".join(cells) + r" \\")
    heads = " & ".join(rf"\multicolumn{{3}}{{c}}{{{COUNTRIES[iso]}}}" for iso in order)
    rules = " ".join(rf"\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}" for i in range(3))
    sub = " & ".join(f"Annual & Monthly & {national_name[iso]}" for iso in order)
    path.write_text(
        "\\begin{table}[htbp]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{4pt}\n"
        f"\\caption{{Reported measles cases, {FIRST_YEAR}--{LAST_YEAR}: WHO annual series (Joint Reporting Form; "
        "the series OWID republishes), sum of WHO's monthly surveillance data (provisional data based on monthly "
        f"data reported to WHO (Geneva) as of {monthly_as_of}), and national totals (CDC, PHAC, DGE). "
        "-- means no value in the annual series, or no monthly records for that year.}\n"
        "\\label{tab:whocheck}\n"
        "\\begin{tabular}{r rrr rrr rrr}\n\\toprule\n"
        f" & {heads} \\\\\n{rules}\nYear & {sub} \\\\\n\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n",
        encoding="utf-8",
        newline="\n",
    )


def write_post2024(annual: pd.DataFrame, monthly: pd.DataFrame, path_csv: Path, path_tex: Path,
                   monthly_as_of: str, data_dir: Path) -> pd.DataFrame:
    """Preliminary totals for the years after the analysis window, one row per country and year."""
    cdc = pd.read_csv(data_dir / CDC_FILE).rename(columns={"Year": "year", "Number of measles cases": "cdc"})
    cdc = cdc.assign(iso3="USA", cdc_note=cdc["Number of measles cases (Annotations)"].fillna(""))
    dge = pd.read_csv(data_dir / NATIONAL_RECENT)[["iso3", "year", "period_end", "cases"]].rename(columns={"cases": "dge"})
    grid = pd.DataFrame([(iso, y) for iso in COUNTRIES for y in RECENT_YEARS], columns=["iso3", "year"])
    df = (grid.merge(annual, how="left").merge(monthly, how="left")
          .merge(cdc[["iso3", "year", "cdc", "cdc_note"]], how="left").merge(dge, how="left"))
    df["national"] = df["cdc"].fillna(df["dge"])
    df["national_source"] = df["iso3"].map({"USA": "CDC (OWID US series, 20261101 snapshot)",
                                            "MEX": "DGE (informe diario, 10 Aug 2026)"}).where(df["national"].notna())
    df["national_period_end"] = df["period_end"].where(df["iso3"] == "MEX")
    # OWID repeats its annotation on every row; keep it only for the year it names
    df["national_note"] = [n if r.iso3 == "USA" and str(r.year) in n else "" for r, n in zip(df.itertuples(), df["cdc_note"].fillna(""))]
    cols = ["iso3", "year", "who_annual", "who_monthly", "first_month", "last_month", "national",
            "national_source", "national_period_end", "national_note"]
    df = df[cols].sort_values(["iso3", "year"]).reset_index(drop=True)
    for c in ["who_annual", "who_monthly", "national"]:
        df[c] = df[c].astype("Int64")
    df.to_csv(path_csv, index=False, lineterminator="\n")

    month = lambda m: pd.Timestamp(f"2000-{m}-01").strftime("%b")

    def num(v) -> str:
        if pd.isna(v):
            return "--"
        return f"{int(v):,}".replace(",", "\\,") if v >= 10000 else f"{int(v)}"

    lines = []
    for _, r in df.iterrows():
        monthly_cell = num(r.who_monthly)
        if pd.notna(r.who_monthly) and (r.first_month, r.last_month) != ("01", "12"):
            monthly_cell += f" ({month(r.first_month)}--{month(r.last_month)})"
        national_cell = num(r.national)
        if r.iso3 == "MEX" and pd.notna(r.national) and r.national_period_end != f"{r.year}-12-31":
            end = pd.Timestamp(r.national_period_end)
            national_cell += f" (to {end.day} {end.strftime('%b')})"
        if r.iso3 == "USA" and str(r.year) in r.national_note:
            national_cell += " (incomplete)"
        lines.append(f"{COUNTRIES[r.iso3]} & {r.year} & {num(r.who_annual)} & {monthly_cell} & {national_cell} \\\\")
    path_tex.write_text(
        "\\begin{table}[htbp]\n\\centering\n\\small\n"
        "\\caption{Reported measles cases after the analysis window, 2025 and part of 2026 (preliminary). "
        "WHO annual series (Joint Reporting Form); sum of WHO's monthly surveillance data (provisional data based on "
        f"monthly data reported to WHO (Geneva) as of {monthly_as_of}); national totals: CDC through the OWID US "
        "series (snapshot of 1 November 2026; 2026 flagged incomplete) and, for Mexico, the DGE outbreak report of "
        "10 August 2026, whose 2026 value is its total since 1 January 2025 minus its 2025 total. No national total "
        "for Canada was recorded; the outbreak that began there in October 2024 exceeded 5000 confirmed cases. "
        "-- means no value.}\n"
        "\\label{tab:post2024}\n"
        "\\begin{tabular}{llrrr}\n\\toprule\n"
        "Country & Year & WHO annual & WHO monthly & National \\\\\n\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n",
        encoding="utf-8",
        newline="\n",
    )
    return df


def monthly_as_of(path: Path) -> str:
    """Release month from the monthly file's Read Me sheet (a 'YYYY-MM' cell), e.g. 'October 2026'."""
    readme = pd.read_excel(path, sheet_name="Read Me", header=None, dtype=str)
    stamps = [v for v in readme.stack() if isinstance(v, str) and len(v) == 7 and v[4] == "-" and v[:4].isdigit()]
    if len(stamps) != 1:
        raise ValueError(f"Expected one YYYY-MM release stamp in {path} Read Me sheet, found {stamps}")
    return pd.Timestamp(stamps[0] + "-01").strftime("%B %Y")


def validate_inputs(snapshot: str, data_dir: Path) -> tuple[Path, Path]:
    """Check every input against its manifest in data_dir; return the monthly and annual WHO files."""
    monthly_p = data_dir / f"who_measles_monthly_{snapshot}.xlsx"
    annual_p = data_dir / f"who_measles_annual_{snapshot}.csv"
    validate(data_dir / f"data_provenance_who_{snapshot}.csv", [monthly_p, annual_p])
    validate(data_dir / f"data_provenance_{CDC_SNAPSHOT}.csv", [data_dir / CDC_FILE], exact=False)
    validate(data_dir / "data_provenance_national.csv", [data_dir / NATIONAL])
    validate(data_dir / "data_provenance_national_2025_2026.csv", [data_dir / NATIONAL_RECENT])
    return monthly_p, annual_p


def main(snapshot: str = "20261107", data_dir: Path = DATA, out_dir: Path = RESULTS) -> None:
    monthly_p, annual_p = validate_inputs(snapshot, data_dir)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)

    annual, monthly = load_annual(annual_p), load_monthly(monthly_p)
    grid = pd.DataFrame(
        [(iso, y) for iso in COUNTRIES for y in range(FIRST_YEAR, LAST_YEAR + 1)], columns=["iso3", "year"]
    )
    df = (
        grid.merge(annual, how="left")
        .merge(monthly, how="left")
        .merge(load_national(data_dir), how="left")
    )
    df = df.sort_values(["iso3", "year"]).reset_index(drop=True)
    prev = df.groupby("iso3")["national"].shift(1)
    df["annual_eq_national_same_year"] = [within_one(a, n) for a, n in zip(df["who_annual"], df["national"])]
    df["annual_eq_national_prev_year"] = [within_one(a, n) for a, n in zip(df["who_annual"], prev)]
    df["annual_eq_monthly_same_year"] = [within_one(a, m) for a, m in zip(df["who_annual"], df["who_monthly"])]
    df["monthly_eq_national_same_year"] = [within_one(m, n) for m, n in zip(df["who_monthly"], df["national"])]

    cols = [
        "iso3", "year", "who_annual", "who_monthly", "months_reported", "national",
        "annual_eq_national_same_year", "annual_eq_national_prev_year",
        "annual_eq_monthly_same_year", "monthly_eq_national_same_year", "national_source",
    ]
    out = out_dir / "who_series_check.csv"
    df[cols].to_csv(out, index=False, lineterminator="\n")
    tex = out_dir / "tables" / "who_check.tex"
    write_tex_table(df, tex, monthly_as_of(monthly_p))
    print(f"Wrote {out} and {tex}\n")
    recent = write_post2024(annual, monthly, out_dir / "post2024_totals.csv", out_dir / "tables" / "post2024.tex",
                            monthly_as_of(monthly_p), data_dir)
    print(recent.to_string(index=False), "\n")

    show = df[cols[:6]].copy()
    for c in ["who_annual", "who_monthly", "months_reported", "national"]:
        show[c] = show[c].astype("Int64")
    print(show.to_string(index=False))
    print("\nYears (of those with a WHO annual value) where WHO annual is within one case of:")
    for iso, g in df[df["who_annual"].notna()].groupby("iso3"):
        print(
            f"  {COUNTRIES[iso]}: national same year {g['annual_eq_national_same_year'].sum()}/{len(g)}, "
            f"national previous year {g['annual_eq_national_prev_year'].sum()}/{len(g)}, "
            f"WHO monthly same year {g['annual_eq_monthly_same_year'].sum()}/{len(g)}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--snapshot", default="20261107", help="Download date of the WHO files (YYYYMMDD)")
    parser.add_argument("--data-dir", type=Path, default=DATA, help="Folder with the input files and manifests")
    parser.add_argument("--out-dir", type=Path, default=RESULTS, help="Folder for the outputs")
    args = parser.parse_args()
    main(args.snapshot, args.data_dir, args.out_dir)
