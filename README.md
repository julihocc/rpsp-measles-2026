# Measles resurgence in Canada, Mexico and the United States, 2025-2026

Code and data to reproduce every number, table and figure of the article "Measles Resurgence in Canada, Mexico and the United States, 2025–2026, from Open Surveillance Data" and of its supplementary material.

The analysis uses archived copies of openly available files: WHO's provisional monthly measles surveillance data, WHO's annual reported cases, Our World in Data (OWID) series of reported cases, vaccine coverage and population, and national totals compiled from the public health agencies of the three countries. Nothing is downloaded when the analysis runs.

## Requirements

[uv](https://docs.astral.sh/uv/). It installs Python 3.12.10 and the exact package versions in `uv.lock`.

## Reproduce

```bash
uv sync
uv run python scripts/analyze_north_america.py
uv run python scripts/check_who_series.py
uv run python scripts/analyze_resurgence.py
```

Run the scripts in this order: `analyze_resurgence.py` checks its 1980s and 2010s rates against the output of `analyze_north_america.py`. Each script first checks every input file against the SHA-256 checksums in `data/data_provenance_*.csv` and stops if one does not match.

The scripts overwrite the files in `results/`. If the run reproduces the published results, `git status --porcelain results/` prints nothing.

## Outputs

| Article | File in `results/` | Script |
|---|---|---|
| Table 1 | `tables/resurgence_rates.tex` | `analyze_resurgence.py` |
| Table 2 | `tables/resurgence_phases.tex` | `analyze_resurgence.py` |
| Table 3 | `tables/resurgence_coverage.tex` | `analyze_resurgence.py` |
| Figure 1 | `figures/fig_cases_per100k.pdf` | `analyze_north_america.py` |
| Figure 2 | `figures/fig_monthly.pdf` | `analyze_resurgence.py` |
| Table S2 | `tables/post2024.tex`, `post2024_totals.csv` | `check_who_series.py` |
| Table S3 | `tables/who_check.tex`, `who_series_check.csv` | `check_who_series.py` |
| Table S4 | `tables/association.tex` | `analyze_north_america.py` |
| Table S5 | `tables/robustness.tex` | `analyze_north_america.py` |
| Section S2 | `results.json` | `analyze_north_america.py` |
| Text of Results | `results_resurgence.json` | `analyze_resurgence.py` |

Table S1 lists the input files in `data/` with their checksums.

## Data

| File in `data/` | Source | Accessed |
|---|---|---|
| `who_measles_monthly_20261107.xlsx` | WHO, provisional monthly measles and rubella data, immunization data portal | 7 November 2026 |
| `who_measles_annual_20261107.csv` | WHO, Global Health Observatory indicator WHS3_62 (reported measles cases) | 7 November 2026 |
| `reported_cases_measles_20261101.csv` | OWID, reported-cases-of-measles (WHO data) | 1 November 2026 |
| `measles_us_cases_20261101.csv` | OWID, number-of-measles-cases (United States, CDC data) | 1 November 2026 |
| `measles_vaccination_mcv1_20261101.csv`, `measles_vaccination_mcv2_20261101.csv` | OWID, measles vaccine coverage (WHO and UNICEF estimates) | 1 November 2026 |
| `population_20261107.csv`, `population_un_projections_20261107.csv` | OWID, population (UN World Population Prospects 2024) | 7 November 2026 |
| `national_measles_totals.csv` | Annual totals, 2015-2024, from the Public Health Agency of Canada (PHAC), Mexico's Dirección General de Epidemiología (DGE) and CDC; one cited source and URL per row | compiled 7 November 2026 |
| `national_measles_totals_2025_2026.csv` | Preliminary 2025-2026 totals for Mexico from the DGE; source and cut-off date per row | compiled 7 November 2026 |

The files are kept byte for byte as downloaded (`.gitattributes` turns off line-ending conversion in `data/`). The manifests `data_provenance_*.csv` give the size and SHA-256 checksum of each file.

Terms of use:

- OWID data are published under CC BY 4.0 (https://ourworldindata.org).
- WHO data are copyright of the World Health Organization and are reused under WHO's terms and conditions of use for WHO data compilations. They were reported to WHO by the national authorities of the countries concerned. Their use implies no endorsement by WHO.
- National totals are public figures from PHAC (https://health-infobase.canada.ca/measles-rubella/), the DGE (https://www.gob.mx/salud) and CDC (through OWID); each row cites its source.

`scripts/download_owid_data.py` and `scripts/download_who_data.py` download new copies of the OWID and WHO files into `data/raw/`, which git ignores. WHO and OWID revise these files without versioning them, so a new download will generally differ from the archived copies.

## License

The code is released under the MIT License (see `LICENSE`). The data keep the terms of their sources, listed above.

## Citation

If you use this code or data, please cite the article. The full reference will be added here on publication.
