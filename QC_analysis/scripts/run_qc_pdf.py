#!/usr/bin/env python
"""Create the selected Figure 2D-F QC panels as three standalone PDFs."""

from __future__ import annotations

import argparse
from pathlib import Path

from qc_pdf_report import create_figure_pdfs
from qc_analysis import (
    build_combined_matrix,
    load_all_reports,
    pearson_corr_matrix,
    prepare_complete_matrix,
    qc_cv_table,
    sample_qc_summary_from_reports,
)

PACKAGE_DIR = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DIA-LFQ QC analysis and export three standalone Figure 2 PDFs.")
    parser.add_argument(
        "--data-dir",
        default=str(PACKAGE_DIR / "data"),
        help="Directory containing the G1-G6 protein-group CSV reports.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(PACKAGE_DIR / "outputs"),
        help="Directory for the three standalone PDF panels.",
    )
    parser.add_argument("--min-detection-fraction", type=float, default=0.8, help="Minimum detection fraction for Figure 2F.")
    args = parser.parse_args()

    reports = load_all_reports(Path(args.data_dir))
    matrix, sample_meta, annotation = build_combined_matrix(reports)
    sample_qc = sample_qc_summary_from_reports(reports)
    qc_cv = qc_cv_table(matrix, sample_meta, annotation)
    complete = prepare_complete_matrix(matrix, args.min_detection_fraction)
    corr = pearson_corr_matrix(complete)
    out_dir = Path(args.out_dir)
    for obsolete_name in ["qc_main_figures.pdf", "figure_2D_2E_2F_qc.pdf"]:
        (out_dir / obsolete_name).unlink(missing_ok=True)
    outputs = create_figure_pdfs(sample_qc, qc_cv, corr, out_dir)

    print(f"Loaded reports: {len(reports)}")
    print(f"Combined matrix: {matrix.shape[0]} proteins x {matrix.shape[1]} samples")
    print(f"Filtered matrix for correlation: {complete.shape[0]} proteins x {complete.shape[1]} samples")
    for output in outputs:
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
