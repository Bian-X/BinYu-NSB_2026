#!/usr/bin/env python
"""QC drift curve for the reviewer response.

Each pooled-QC injection's Pearson correlation to the across-QC median
profile, plotted against its chronological injection order, with a
descriptive linear fit. No confidence band or inferential P value is drawn
because all observations share the across-QC median reference.

Data preprocessing is IDENTICAL to Figure 2F (manuscript Section 3.10.3):
  combined protein-group matrix -> >=80% detection filter -> log2 ->
  sample-median normalization -> protein-median imputation.
The correlation is then computed on the 25 QC columns only.

Outputs (written to --out-dir):
  figure_QC_drift.pdf / .png   publication-style figure
  qc_drift_metrics.csv         per-QC correlation + injection order

Usage:
  python qc_drift_curve.py --data-dir ../data \
      --metadata ../Sample_metadata.xlsx --out-dir outputs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# The module sits next to this script when installed in
# QC_analysis/scripts; the fallback keeps the development workspace
# layout working.
_SCRIPT_DIR = Path(__file__).resolve().parent
for _candidate in (
    _SCRIPT_DIR,
    _SCRIPT_DIR.parent / "BinYu-NSB_2026" / "QC_analysis" / "scripts",
):
    if (_candidate / "qc_analysis.py").exists():
        sys.path.insert(0, str(_candidate))
        break
else:
    raise FileNotFoundError(
        "qc_analysis.py not found; place this script in QC_analysis/scripts "
        "or keep the BinYu-NSB_2026 workspace layout"
    )

import qc_analysis as qa  # noqa: E402
from metadata_io import read_sample_metadata, match_metadata


def build_qc_drift_table(data_dir: Path, metadata_path: Path) -> pd.DataFrame:
    reports = qa.load_all_reports(data_dir)
    matrix, _sample_meta, _annotation = qa.build_combined_matrix(reports)

    imputed = qa.prepare_complete_matrix(matrix, min_detection_fraction=0.8)
    print(f"Matrix after Figure 2F preprocessing: {imputed.shape[0]} protein groups x {imputed.shape[1]} samples")

    meta = read_sample_metadata(metadata_path)
    key = match_metadata(meta, imputed.columns)
    qc_meta = key.loc[key["sample_type"] == "QC"]

    qc_cols = [c for c in imputed.columns if c in qc_meta.index]
    if len(qc_cols) < 3:
        raise ValueError("At least three matched QC injections are required")
    qc_mat = imputed[qc_cols]
    qc_median = qc_mat.median(axis=1)
    corr = qc_mat.corrwith(qc_median)  # Pearson r of each QC vs QC median

    table = pd.DataFrame({
        "sample_id": qc_cols,
        "injection_order": qc_meta.loc[qc_cols, "injection_order"].astype(int).values,
        "correlation_to_qc_median": corr.values,
    }).sort_values("injection_order").reset_index(drop=True)
    return table


def plot_drift(table: pd.DataFrame, out_base: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    })

    x = table["injection_order"].values.astype(float)
    y = table["correlation_to_qc_median"].values

    fit = stats.linregress(x, y)
    x_line = np.array([x.min(), x.max()])
    y_line = fit.intercept + fit.slope * x_line

    n = len(x)

    fig, ax = plt.subplots(figsize=(4.6, 3.0))

    ax.plot(x_line, y_line, color="#08519C", lw=1.2, zorder=2,
            label="Linear fit")
    ax.scatter(x, y, s=16, facecolor="#08519C", edgecolor="white",
               linewidth=0.5, zorder=3, label=f"Pooled QC (n = {n})")

    order_span = x.max() - x.min()
    drift_total = fit.slope * order_span
    stats_text = (
        f"r range: {y.min():.4f}\u2013{y.max():.4f}\n"
        f"Slope: {fit.slope:+.2e} r / injection\n"
        f"Correlation change over {order_span:.0f} order units: {drift_total:+.4f}\n"
        f"Descriptive $R^2$ = {fit.rvalue ** 2:.4f}; shared reference"
    )
    ax.text(0.03, 0.04, stats_text, transform=ax.transAxes, fontsize=6.5,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#999999", linewidth=0.5, alpha=0.9))

    ax.set_xlabel("LC-MS/MS injection order", fontsize=8)
    ax.set_ylabel("Pearson r vs QC median", fontsize=8)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    bounds = np.concatenate([y, y_line])
    pad = max(0.001, 0.08 * np.ptp(bounds))
    ax.set_ylim(bounds.min() - pad, bounds.max() + pad)
    ax.tick_params(direction="out", length=2.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.17), ncol=3,
              fontsize=6.5, frameon=False, borderaxespad=0.2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.savefig(out_base.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_base.with_suffix('.pdf')} and {out_base.with_suffix('.png')}")

    print(f"slope = {fit.slope:.3e} r per order unit; correlation change = {drift_total:+.4f}")
    print("Correlation to a shared QC reference is descriptive and does not bound absolute signal drift.")
    print(f"Descriptive R^2 = {fit.rvalue ** 2:.4f}; no inferential regression P is reported")
    print(f"r range: {y.min():.4f} - {y.max():.4f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, required=True, help="Folder with the 12 Spectronaut CSV reports")
    parser.add_argument("--metadata", type=Path, required=True, help="Sample_metadata.xlsx with injection_order")
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    table = build_qc_drift_table(args.data_dir, args.metadata)
    table.to_csv(args.out_dir / "qc_drift_metrics.csv", index=False)
    plot_drift(table, args.out_dir / "figure_QC_drift")


if __name__ == "__main__":
    main()
