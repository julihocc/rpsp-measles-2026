"""
Download WHO measles surveillance data used to check the WHO-derived annual series.

Two files are saved to data/raw/who/ (or --out-dir) with the download date (YYYYMMDD) in the name:

- who_measles_monthly_<date>.xlsx: WHO's provisional monthly measles and rubella
  surveillance file (cases by country, year, month and final classification since
  2012), from the WHO immunization data portal. Monthly data exist only for
  countries that report monthly or weekly.
- who_measles_annual_<date>.csv: WHO's official annual reported measles cases
  (Global Health Observatory indicator WHS3_62, from the WHO/UNICEF Joint Reporting
  Form), all countries and years. This is the series that Our World in Data
  republishes as reported-cases-of-measles, and it matches the immunization data
  portal's reported-cases export.

WHO revises both sources and does not version them, so a new download cannot
reproduce an earlier snapshot; the copies used in the analysis are archived in data/.
"""

import argparse
import csv
from datetime import datetime
from pathlib import Path

import requests

MONTHLY_URL = (
    "https://immunizationdata.who.int/docs/librariesprovider21/measles-and-rubella/"
    "404-table-web-epi-curve-data.xlsx"
)
GHO_URL = "https://ghoapi.azureedge.net/api/WHS3_62"
HEADERS = {"User-Agent": "measles-research-download-script"}
ROOT = Path(__file__).resolve().parents[1]
GHO_COLUMNS = ["iso3", "year", "cases", "spatial_type", "last_updated"]


def download_monthly(output_path: Path) -> None:
    """Save WHO's provisional monthly measles and rubella file unchanged."""
    response = requests.get(MONTHLY_URL, headers=HEADERS, timeout=300)
    response.raise_for_status()
    output_path.write_bytes(response.content)


def download_annual(output_path: Path) -> None:
    """Save the GHO annual reported measles cases (all countries) as CSV."""
    response = requests.get(GHO_URL, headers=HEADERS, timeout=300)
    response.raise_for_status()
    rows = [
        {
            "iso3": r["SpatialDim"],
            "year": r["TimeDim"],
            # WHO lists some country-years with no value; keep them as empty cells
            "cases": "" if r["NumericValue"] is None else int(r["NumericValue"]),
            "spatial_type": r["SpatialDimType"],
            "last_updated": (r.get("Date") or "")[:10],
        }
        for r in response.json()["value"]
    ]
    rows.sort(key=lambda r: (r["spatial_type"], r["iso3"], r["year"]))
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GHO_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(output_dir: Path = ROOT / "data" / "raw" / "who") -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d")

    downloads = {
        f"who_measles_monthly_{timestamp}.xlsx": download_monthly,
        f"who_measles_annual_{timestamp}.csv": download_annual,
    }
    failed = []
    for name, download in downloads.items():
        print(f"Downloading {name}...")
        try:
            download(output_dir / name)
            print(f"  Saved to {output_dir / name}")
        except requests.exceptions.RequestException as e:
            print(f"  Error downloading {name}: {e}")
            failed.append(name)

    if failed:
        raise SystemExit(f"Failed downloads: {', '.join(failed)}")
    print("\nSource: World Health Organization, immunization data portal and Global Health Observatory.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data" / "raw" / "who",
                        help="Folder for the downloaded files")
    main(parser.parse_args().out_dir)
