# Figure 3: G4 versus G2 proteomics workflow

Reproducible scripts for the five Figure 3 data-usage demonstrations from a
DIA label-free quantitative (DIA-LFQ) mouse brain proteomics dataset. The
workflow compares paired tone-shock learning (G4) with temporally unpaired
tone-shock exposure (G2) across 11 brain regions.

## What this repository contains

- R scripts for differential abundance, Gene Ontology enrichment, heatmaps,
  and regional log2-fold-change correlations.
- A Python script for cross-region principal component analysis (PCA).
- Windows PowerShell and Linux/macOS Bash runners.
- Dependency lists and a recorded R session information file.
- GitHub Actions checks that validate script syntax without requiring data.

The repository deliberately **does not contain raw mass-spectrometry files,
processed input tables, sample metadata, or generated figures/results**. These
files can be large and may be subject to data-release restrictions.

## Required input data

Place the 11 non-QC Spectronaut-style region CSV reports in `data/`, or pass
their folder to a runner. File names must contain each of these region codes:

`PMD`, `LPB`, `VHPC`, `PVH`, `BLA`, `DHPC`, `CEA`, `MSC`, `SSC`, `IL`, `PL`.

The R analyses require exactly four G2 and four G4 quantity columns per
region. The PCA script expects 24 quantity columns (G1-G6, four biological
replicates per group). See [data/README.md](data/README.md) for the expected
columns and naming pattern.

## Environment setup

Tested environment: Python 3.12.13; R 4.4.3 with Bioconductor 3.20. The
direct Python dependencies are pinned in `environment/requirements.txt`; the
complete tested Python dependency set is recorded in
`environment/python-requirements-lock.txt`.

```bash
python -m pip install -r environment/requirements.txt
Rscript environment/install_R_packages.R
```

The R package list is in `environment/r-packages.txt`; exact direct-package
versions are recorded in `environment/r-package-versions.csv`, and the full R
execution environment is in `environment/sessionInfo.txt`. To use the complete
Python lock file instead of the direct dependencies, run
`python -m pip install -r environment/python-requirements-lock.txt`.

## Run the complete workflow

Windows PowerShell:

```powershell
.\run_all.ps1 -DataDir "C:\path\to\region_csvs"
```

Linux/macOS:

```bash
DATA_DIR=/path/to/region_csvs bash run_all.sh
```

The runners write PDF figures to `figures/` and result tables to `results/`.

## Generated figures

| File | Description |
| --- | --- |
| `BLA_volcano_G4_vs_G2.pdf` | BLA G4/G2 volcano plot |
| `BLA_GO_enrichment_G4_vs_G2.pdf` | BLA GO enrichment bubble plot |
| `brain_region_pca.pdf` | PCA of samples from 11 brain regions |
| `multi_region_G2_vs_G4_heatmap.pdf` | Region-balanced G2/G4 protein heatmap |
| `brain_region_log2FC_correlation.pdf` | Correlation heatmap of regional G4-G2 log2FC profiles |

## Statistical interpretation

The differential-protein and GO demonstrations use raw P-value thresholds and
are explicitly exploratory. They must not be interpreted as confirmatory
findings without multiple-testing adjustment and an analysis plan appropriate
to the intended claim.
