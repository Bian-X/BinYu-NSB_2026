#!/usr/bin/env python
"""Figure 2F redrawn in LC-MS/MS injection order.

This script reproduces the manuscript's Figure 2F sample-correlation
methodology exactly (see BinYu-NSB_2026/QC_analysis and manuscript
Section 3.10.3), but orders the 289 samples left-to-right by the
chronological LC-MS/MS injection order recorded in
Sample_metadata.xlsx (field `injection_order`) instead of the
region-blocked order used in the published panel.

Figure 2F methodology (unchanged from the manuscript):
  1. Combine the 12 Spectronaut protein-group reports into one
     protein-group x sample matrix (11,322 x 289).
  2. Retain protein groups detected in at least 80% of all 289
     measurements (yields 8,072 protein groups).
  3. log2-transform the quantities.
  4. Sample-median normalize (subtract each sample's median).
  5. Impute remaining missing values with the protein-wise median.
  6. Compute the pairwise Pearson correlation matrix across samples.

Only the sample ordering differs: columns/rows are sorted by
`injection_order` (1 = first injection). Pooled-QC injections are
interspersed through the acquisition series (2 at the start, usually 1 after
12 biological samples, 2 at the end). In the retained sequence there are 8
biological samples between QC_12 and QC_13 and 16 between QC_13 and QC_14;
all positions are read from metadata. This view allows inspection of temporal
patterns that the region-blocked view can obscure.

Outputs (written to --out-dir):
  figure_2F_injection_order.pdf / .png  heatmap
  injection_order_correlation_matrix.csv  the 289 x 289 matrix
  injection_order_sample_key.csv  sample -> injection order mapping

Usage:
  python figure_2F_injection_order.py \
      --data-dir ../data \
      --metadata ../Sample_metadata.xlsx \
      --out-dir outputs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the repository's verified loading / normalization code so the
# matrix construction is identical to the published workflow. The module
# sits next to this script when installed in QC_analysis/scripts; the
# fallback keeps the development workspace layout working.
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

# Region colors reused from the published figure for visual consistency.
REGION_COLORS = qa.REGION_COLORS


def build_figure_2f_correlation(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reproduce the Figure 2F correlation matrix (region order)."""
    reports = qa.load_all_reports(data_dir)
    matrix, sample_meta, _annotation = qa.build_combined_matrix(reports)

    # Every QC view uses exactly one shared preprocessing implementation.
    imputed = qa.prepare_complete_matrix(matrix, min_detection_fraction=0.8)
    n_kept = len(imputed)
    corr = qa.pearson_corr_matrix(imputed)
    print(f"Combined matrix: {matrix.shape[0]} protein groups x {matrix.shape[1]} samples")
    print(f"Retained for correlation (>=80% detection): {n_kept} protein groups")
    return corr, sample_meta


def attach_injection_order(corr: pd.DataFrame, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sort the correlation matrix by chronological injection order."""
    meta = read_sample_metadata(metadata_path)
    key = meta.set_index("sample_id")[["injection_order", "region", "group", "sample_type"]]

    missing = [c for c in corr.columns if c not in key.index]
    if missing:
        raise ValueError(f"{len(missing)} matrix samples lack metadata rows, e.g. {missing[:3]}")

    key = match_metadata(meta, corr.columns)
    order = key.sort_values("injection_order").index
    return corr.loc[order, order], key.loc[order]


def plot_heatmap(corr: pd.DataFrame, key: pd.DataFrame, out_base: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch

    n = corr.shape[0]
    regions = key["region"].fillna("QC")
    colors = [REGION_COLORS.get(r, "#000000") for r in regions]

    cmap = LinearSegmentedColormap.from_list("corr", ["#f7fbff", "#08306b"])
    fig_size = 11
    fig = plt.figure(figsize=(fig_size + 3.2, fig_size + 1.2))
    # [left, bottom, width, height]
    ax_top = fig.add_axes([0.06, 0.785, 0.72, 0.015])  # top color band
    ax_left = fig.add_axes([0.04, 0.06, 0.015, 0.72])  # left color band
    ax = fig.add_axes([0.06, 0.06, 0.72, 0.72])
    cax = fig.add_axes([0.82, 0.35, 0.02, 0.25])

    im = ax.imshow(corr.values, cmap=cmap, vmin=0.85, vmax=1.0, aspect="auto")

    # Region color bands.
    band = np.arange(n)
    ax_top.scatter(band, np.zeros(n), c=colors, marker="s", s=4, linewidths=0)
    ax_top.set_xlim(-0.5, n - 0.5)
    ax_top.axis("off")
    ax_left.scatter(np.zeros(n), band, c=colors, marker="s", s=4, linewidths=0)
    ax_left.set_ylim(n - 0.5, -0.5)
    ax_left.axis("off")

    # QC injection positions: small tick marks just under the heatmap
    # (no lines across the heatmap body).
    qc_pos = np.where(key["sample_type"].values == "QC")[0]
    tick_len = n * 0.012
    for pos in qc_pos:
        ax.plot([pos, pos], [n - 0.5, n - 0.5 + tick_len],
                color="black", lw=0.8, solid_capstyle="butt", clip_on=False)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Samples sorted by LC-MS/MS injection order (1 → 289); black ticks below = pooled-QC injections")
    fig.suptitle(
        "Figure 2F reordered by injection order\n"
        "Pearson correlation of 8,072 protein groups (log2, sample-median normalized, protein-median imputed)",
        fontsize=11,
        y=0.98,
    )

    fig.colorbar(im, cax=cax, label="Pearson r")
    present = [r for r in dict.fromkeys(regions) ]
    handles = [Patch(facecolor=REGION_COLORS.get(r, "#000000"), label=r) for r in REGION_COLORS if r in set(present)]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.92, 0.72), frameon=False, title="Region")

    fig.savefig(out_base.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_base.with_suffix('.pdf')} and {out_base.with_suffix('.png')}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, required=True, help="Folder with the 12 Spectronaut CSV reports")
    parser.add_argument("--metadata", type=Path, required=True, help="Sample_metadata.xlsx with injection_order")
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    corr, _sample_meta = build_figure_2f_correlation(args.data_dir)
    corr_ordered, key = attach_injection_order(corr, args.metadata)

    corr_ordered.to_csv(args.out_dir / "injection_order_correlation_matrix.csv")
    key.to_csv(args.out_dir / "injection_order_sample_key.csv")
    plot_heatmap(corr_ordered, key, args.out_dir / "figure_2F_injection_order")

    # Quick drift check: correlation of each sample with the first vs last QC.
    qc_ids = key.index[key["sample_type"] == "QC"]
    first_qc, last_qc = qc_ids[0], qc_ids[-1]
    print(f"First QC injection: {first_qc} (order {key.loc[first_qc, 'injection_order']})")
    print(f"Last QC injection:  {last_qc} (order {key.loc[last_qc, 'injection_order']})")
    print(f"Correlation first-vs-last QC: {corr_ordered.loc[first_qc, last_qc]:.4f}")


if __name__ == "__main__":
    main()
