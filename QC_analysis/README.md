# DIA-LFQ proteomics QC figures

Reproducible code for the three proteomics quality-control panels in Figure 2:

- **Figure 2D**: protein-group detection depth and within-report missingness across 289 LC-MS/MS measurements.
- **Figure 2E**: pooled-QC protein-group coefficient of variation (CV) versus mean log2 abundance.
- **Figure 2F**: Pearson sample-correlation heatmap.

The repository contains code and documentation only. The 12 original DIA-LFQ protein-group reports, derived tables, and publication figures are deliberately excluded from version control.

## Repository layout

```text
scripts/       analysis module and PDF entry point
data/          place the 12 source CSV reports here (not tracked)
outputs/       generated PDFs (not tracked)
environment/   verified Python and dependency versions
docs/          methods and manuscript-alignment record
```

## Requirements

The release was verified with CPython 3.12.13 on Windows and the pinned packages in `requirements.txt`:

```text
numpy==2.3.5
pandas==3.0.1
reportlab==4.4.9
```

Create an isolated environment and install the requirements:

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS/Linux
# source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Input data

Copy the 12 original protein-group CSV reports to `data/`. There must be one pooled-QC report and eleven regional reports. Each report must include `PG.ProteinAccessions` and one or more columns ending in `.PG.Quantity`; sample names must follow the original `Astral_MSP2600441_...raw` convention. See [data/README.md](data/README.md) for details.

The raw/processed proteomics reports are not included because their sharing and accession conditions must be determined by the study authors. Do not upload them to a public repository unless their release is authorized.

## Run

From the repository root:

```bash
python scripts/run_qc_pdf.py --data-dir data --out-dir outputs
```

The command writes these three one-page PDFs:

```text
outputs/figure_2D_detection_missing.pdf
outputs/figure_2E_qc_cv_abundance.pdf
outputs/figure_2F_sample_correlation.pdf
```

For the archived input reports, the expected console summary is 12 reports, 11,322 protein groups x 289 samples, and 8,072 protein groups x 289 samples after the 80% detection filter.

## Methods

- **Figure 2D:** counts and missingness are calculated within each report's own protein-group universe, without imputation.
- **Figure 2E:** proteins quantified in at least 20 of 25 pooled-QC injections are retained. CV is `100 x sample SD / mean` on non-log intensities; missing values are not imputed. Values above 100% are displayed at 100% in the plot.
- **Figure 2F:** protein groups detected in at least 80% of all measurements are retained; quantities are log2 transformed, sample-median normalized, then remaining values are imputed with the protein-wise median before Pearson correlation. The heatmap is descriptive because median imputation can inflate correlations.

Detailed correspondence with the revised manuscript is recorded in [docs/MANUSCRIPT_ALIGNMENT.md](docs/MANUSCRIPT_ALIGNMENT.md).

## Reproducibility and scope

`scripts/run_qc_pdf.py` is the supported public entry point. `scripts/qc_analysis.py` also retains the archival SVG routine for historical reproducibility, but it is not the supported release workflow. `environment/python-requirements-lock.txt` records the complete verified Python package set.

No license is included. Add a license only after the code owners choose the intended reuse terms.

