# Multi-region mouse DIA-LFQ proteomics

Analysis code for Figure 2D-F and Figure 3A-E of *A within-animal matched
multi-region proteomics dataset of fear conditioning and memory retrieval in
mice*. Version: **2026.10.26**.

The workflow calculates acquisition quality control, differential protein
abundance, Gene Ontology (GO) enrichment, principal component analysis (PCA),
an exploratory heatmap, and correlations between regional abundance changes.
Each of the eight manuscript panels is exported as PDF, SVG, PNG, and TIFF.
See [METHODS.md](METHODS.md) for preprocessing, statistical models, and
interpretation limits.

## Input data

The analysis requires 12 processed protein-group CSV reports: one for each of
11 brain regions and one pooled-QC report. Obtain the study data from the
[dataset record](https://doi.org/10.57760/sciencedb.43331) and place the reports
in a single directory. Report filenames, byte sizes, and SHA-256 checksums are
listed in [metadata/input_manifest.csv](metadata/input_manifest.csv). The runner
verifies these files before analysis. The QC reader selects these 12 reports
from the manifest, so an unrelated metadata CSV is not treated as a report.

Biological quantity columns must use the original **G1-G6** identifiers, with
four replicates per group and region. For example:

```text
[1] Astral_MSP2600441_1_G1_1_PL.raw.PG.Quantity
```

The annotation fields are `PG.ProteinAccessions`, `PG.Genes`, `PG.Organisms`,
`PG.ProteinDescriptions`, and
`PG.NrOfStrippedSequencesIdentified (Experiment-wide)`. Sample identifiers are
read directly without group relabeling. A space immediately before a report's
`.csv` extension is accepted when its contents match the input manifest.

[metadata/Sample_metadata.xlsx](metadata/Sample_metadata.xlsx) is included and
contains 289 acquisitions: 264 regional biological measurements and 25 pooled-QC
acquisitions. Required fields are `raw_file`, `sample_type`, `region`, `group`,
and `injection_order`. Headers may occupy row 1 or row 2. The reader preserves
missing fields and validates sample identities and injection order; blank QC
region fields are interpreted as `QC`. The workbook is read without modification.
Behavioral columns are outside this analysis; freezing times are in seconds.
An alternative workbook or flat CSV with the required fields can be supplied
with `--metadata`.

## Software requirements

The reference environment uses CPython 3.12.14, R 4.6.1, and Bioconductor 3.23
on Windows. Python 3.12 patch versions are accepted and recorded by the root
runner. Full runs on Linux and macOS have not been validated.

| Python package | Version | Use |
| --- | --- | --- |
| NumPy | 2.5.4 | Numerical arrays |
| pandas | 3.0.1 | Report and metadata tables |
| SciPy | 1.18.1 | Statistics and clustering |
| scikit-learn | 1.9.1 | PCA |
| Matplotlib | 3.11.2 | Plots and panel export |
| openpyxl | 3.1.5 | Excel metadata input |
| ReportLab | 4.4.9 | Supplementary QC PDF report |
| Pillow | 12.3.0 | Raster image export |

| R package | Version | Use |
| --- | --- | --- |
| limma | 3.68.5 | Differential abundance models |
| AnnotationDbi | 1.74.0 | Gene and GO annotation queries |
| org.Mm.eg.db | 3.23.0 | Mouse gene annotations |
| GO.db | 3.23.1 | GO term annotations |
| statmod | 1.5.2 | Robust empirical Bayes moderation |
| Matrix | 1.7.5 | Optional exhaustive screen and complete-case analysis |

Python dependencies are pinned in [requirements.txt](requirements.txt).
The R installer checks the versions above and installs missing or mismatched
packages. Matching archived packages may be needed if repository versions change.
The analysis itself does not download reports or annotations.

## Installation and execution

Use Python 3.12 to create the environment. On Windows, run from the repository
root in PowerShell, with `Rscript` available on `PATH`:

```powershell
python -m venv .venv
$python = Join-Path $PWD '.venv/Scripts/python.exe'
& $python -m pip install -r requirements.txt
$env:R_LIBS_USER = Join-Path $PWD '.R_library'
New-Item -ItemType Directory -Force $env:R_LIBS_USER | Out-Null
$env:LC_ALL = 'English_United States.utf8'
Rscript sample_analysis/environment/install_R_packages.R
./run_all.ps1 -DataDir 'D:/proteomics/data' -Python $python -Rscript Rscript -RLibrary $env:R_LIBS_USER -OutputRoot 'D:/proteomics/results'
```

Replace the data and output paths with your own locations. The equivalent
Python entry point is:

```text
python run_all.py --data-dir /path/to/reports --output-root /path/to/results --rscript Rscript
```

Use `--r-library` if R packages are installed outside the default library paths.
`Rscript sample_analysis/environment/install_R_packages.R --check-only` checks
the active libraries without installation.

Figure 2 can run without R. Select either analysis branch with:

```text
python run_all.py --data-dir /path/to/reports --figures figure2
python run_all.py --data-dir /path/to/reports --figures figure3 --rscript Rscript
```

In PowerShell, the corresponding option is `-Figures figure2` or `-Figures
figure3`. The Bash wrapper accepts `DATA_DIR`, `PYTHON`, `RSCRIPT`, `R_LIBS_USER`,
`METADATA`, `OUTPUT_ROOT`, and `FIGURES`; additional command-line options are
forwarded to the Python runner.

Other options:

- `--selection-sensitivity`: run all 165 region/contrast combinations and the
  BLA complete-case sensitivity analysis. Optional results are checked and
  recorded for the current run only when this option is selected.
- `--font-family`: select a font for regional R/PCA plots and manuscript panels;
  Arial was used for the validated Windows figures. Supplementary QC and GO
  plots use their own defaults.
- `--events`: supply measured acquisition events as CSV with `injection_order`,
  `event_type`, and `description`.
- `--rt-metrics`: supply measured QC retention-time metrics as CSV with
  `sample_id` and `rt_metric`; `--rt-label` sets the axis label.
- `--environment compatible`: allow software-version differences while recording
  them. Result differences are reported; this mode does not enforce agreement
  with manuscript numbers. Input-report checksums remain required.

Without `--output-root`, outputs are written under the repository root.
Output directories are created automatically. Missing event or retention-time
records are not estimated.

Use a fresh output directory for each analysis when possible. In a reused
directory, the run manifest lists only files created or rewritten by the active
branches; unrelated older files remain on disk and are not current-run outputs.
Omitting `--events` removes the event table previously exported by this script.

## Figure map

| Panel | Analysis | Main script |
| --- | --- | --- |
| 2D | Observed detections and missingness, all 289 acquisitions | `QC_analysis/scripts/qc_analysis.py` |
| 2E | Original-scale CV, at least 20 of 25 QC detections | `QC_analysis/scripts/run_qc_pdf.py` |
| 2F | Pearson correlation of the normalized and completed acquisition matrix | `QC_analysis/scripts/figure_2F_injection_order.py` |
| 3A | PL G5 minus G1, six-group limma model | `sample_analysis/scripts/R/volcano_PL_G5_vs_G1.R` |
| 3B | BLA G4 minus G1 GO analysis, unique Entrez mapping | `sample_analysis/scripts/R/go_enrichment_BLA_G4_vs_G1.R` |
| 3C | Centered, unscaled regional PCA | `sample_analysis/scripts/python/brain_region_pca.py` |
| 3D | Exploratory G4 minus G2 model and heatmap | `sample_analysis/scripts/R/region_heatmap_G4_vs_G2.R` |
| 3E | Correlation of regional G5 minus G1 effects | `sample_analysis/scripts/R/region_log2fc_correlation_G5_G1.R` |

`sample_analysis/scripts/python/assemble_manuscript_figures.py` renders all eight
final panels. Figure 3B uses a separate contrast from Figure 3A; its GO terms
do not annotate the PL volcano candidate set.

## Outputs and validation

| Directory | Contents |
| --- | --- |
| `manuscript_figures/` | Eight individual panels in PDF, SVG, PNG, and LZW TIFF |
| `sample_analysis/results/` | Differential abundance, GO, PCA, heatmap, and regional correlation tables |
| `sample_analysis/figures/` | Individual analysis plots and sensitivity plots |
| `QC_analysis/outputs/` | Observed QC metrics, correlation matrices, injection-order summaries, and selection audits |
| `validation/` | Run provenance, exported R matrices, and numerical checks |

Raster manuscript exports use 400 dpi, except Figure 2F, which uses 1200 dpi.
Dense data layers are rasterized in PDF/SVG; labels and axes remain vector.
Figure 2A-C and combined manuscript layouts are outside the workflow's scope.

Independent checks compare all 11
regional R/Python matrices and verify BH adjustment, GO probabilities, PCA,
heatmap completion, and QC correlations. The optional screen adds five
manuscript-result checks. Figure bytes can vary with fonts, rendering libraries,
and PDF timestamps.

Run checks that do not require the processed reports with:

```text
python validation/check_source_manifest.py
python validation/test_regressions.py
Rscript validation/test_regressions.R
```

The runner returns a nonzero exit code on failure. Inspect
`validation/run_manifest.json`, `validation/manuscript_result_checks.json`, and,
for a full run, `validation/numeric_verification.json`. For runs launched through
`run_all.py`, the run manifest records software versions, completed subprocess
commands and exit codes, the 12 processed-report hashes, and execution status.
The Figure 2 branch also records the metadata workbook hash. Output hashes are
collected after successful analysis and checks. Optional event and retention-time
input files are not hashed by this runner. Standalone scripts do not write this
root-workflow manifest.
[SOURCE_MANIFEST.sha256](SOURCE_MANIFEST.sha256) identifies the distributed
source files and metadata workbook. GitHub Actions checks source integrity,
syntax, and regression tests without the processed reports.

## Citation and license

Zimu Zhang, Liangmingmin Zhou, Jingtong Xue, Yiheng Xu, Zixu Zhang,
Shiqian Xiang, Wenkui Lu, Yujie Li, Xiao-Ming Li, and Bin Yu.
*A within-animal matched multi-region proteomics dataset of fear conditioning
and memory retrieval in mice.*

Cite the study and the repository version or commit used for the analysis.
Software citation metadata is provided in [CITATION.cff](CITATION.cff).
The dataset identifier is [10.57760/sciencedb.43331](https://doi.org/10.57760/sciencedb.43331).

The code is distributed under the [BSD 3-Clause License](LICENSE).
Dependencies retain their own licenses. The repository includes sample metadata
and excludes raw acquisitions, processed report CSVs, generated results, and
installed software environments. Data access and reuse conditions are specified
in the dataset record.
