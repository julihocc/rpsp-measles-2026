"""
Download measles-related datasets from Our World in Data.

This script downloads CSV data from Our World in Data's public API
and saves it to data/raw/owid/ (or --out-dir), byte for byte as served.
Pass --only to download a subset, e.g. --only population. The files used in the
analysis are archived in data/ with their checksums.
"""

import argparse
from pathlib import Path

import requests
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]

# Our World in Data measles datasets
OWID_DATASETS = {
    "reported_cases_measles": "https://ourworldindata.org/grapher/reported-cases-of-measles.csv",
    "measles_deaths": "https://ourworldindata.org/grapher/deaths-due-to-measles-gbd.csv",
    "measles_vaccination_mcv1": "https://ourworldindata.org/grapher/share-of-children-vaccinated-against-measles.csv",
    "measles_vaccination_mcv2": "https://ourworldindata.org/grapher/share-of-children-vaccinated-with-mcv2.csv",
    "measles_us_cases": "https://ourworldindata.org/grapher/number-of-measles-cases.csv",
    "measles_us_cases_deaths": "https://ourworldindata.org/grapher/measles-cases-and-death.csv",
    "population": "https://ourworldindata.org/grapher/population.csv",
    "population_un_projections": "https://ourworldindata.org/grapher/population-with-un-projections.csv",
}


def download_dataset(url: str, output_path: Path) -> None:
    """
    Download a dataset from URL and save to output path.

    Args:
        url: URL of the CSV dataset
        output_path: Path where the file should be saved
    """
    print(f"Downloading {output_path.name}...")

    response = requests.get(url)
    response.raise_for_status()

    output_path.write_bytes(response.content)
    print(f"  Saved to {output_path}")


def main(only: list[str] | None = None, output_dir: Path = ROOT / "data" / "raw" / "owid"):
    """Download the Our World in Data datasets (all, or only the named ones)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Add timestamp to track when data was downloaded
    timestamp = datetime.now().strftime("%Y%m%d")

    print("Downloading Our World in Data measles datasets...")
    print(f"Output directory: {output_dir}\n")

    for name, url in OWID_DATASETS.items():
        if only and name not in only:
            continue
        output_file = output_dir / f"{name}_{timestamp}.csv"

        try:
            download_dataset(url, output_file)
        except requests.exceptions.RequestException as e:
            print(f"  Error downloading {name}: {e}")

    print("\nDownload complete!")
    print(f"\nDatasets are licensed under Creative Commons BY 4.0")
    print("Citation: Our World in Data (2025) - https://ourworldindata.org")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", choices=sorted(OWID_DATASETS),
                        help="Download only these datasets")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data" / "raw" / "owid",
                        help="Folder for the downloaded files")
    args = parser.parse_args()
    main(args.only, args.out_dir)
