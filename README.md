# Reproducible proteomics analysis workflows

This repository packages two independent, reproducible analysis workflows from the same study. They are intentionally kept separate: they have different aims, input-data requirements, software dependencies, and output figures.

| Folder | Scope | Main entry point |
|---|---|---|
| [sample_analysis](sample_analysis/) | Regional G2-versus-G4 differential-proteomics analyses and Figure 3 panels | run_all.ps1 (Windows) or run_all.sh (macOS/Linux) |
| [QC_analysis](QC_analysis/) | Figure 2D-F proteomics quality-control panels | python scripts/run_qc_pdf.py --data-dir data --out-dir outputs |

## Important

- Read the README.md inside the relevant subfolder before running its workflow.
- The two subfolders are **not merged** and share no code files or configuration. Install each subfolder's dependencies separately.
- Original proteomics reports, derived result tables, and publication figures are excluded from this public code package. Follow the data instructions in each subfolder; do not upload data without authorization.
- Each subfolder retains its own documentation and dependency lock files. The root GitHub Actions workflow performs syntax and entry-point checks for both workflows without requiring protected input data.

## Upload

This is a single Git repository. From this root directory:

```bash
git add .
git commit -m "Release Figure 2 QC and Figure 3 proteomics workflows"
git remote add origin https://github.com/USERNAME/REPOSITORY.git
git push -u origin main
```
