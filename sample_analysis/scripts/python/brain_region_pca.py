#!/usr/bin/env python3
"""PCA of DIA-LFQ protein-group quantities across brain regions."""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
ANNOTATION_COLUMNS = [
    "PG.ProteinAccessions",
    "PG.Genes",
    "PG.Organisms",
    "PG.ProteinDescriptions",
    "PG.NrOfStrippedSequencesIdentified (Experiment-wide)",
]
REGION_CODES = ("PMD", "LPB", "VHPC", "PVH", "BLA", "DHPC", "CEA", "MSC", "SSC", "IL", "PL")
DEFAULT_EXCLUDED_SAMPLE = "Astral_MSP2600441_253_G6_1_PVH"
REGION_PALETTE = {
    "PMD": "#0072B2",
    "LPB": "#D55E00",
    "VHPC": "#009E73",
    "PVH": "#CC79A7",
    "BLA": "#D62728",
    "DHPC": "#56B4E9",
    "CEA": "#7F3C8D",
    "MSC": "#11A579",
    "SSC": "#3969AC",
    "IL": "#F0E442",
    "PL": "#E73F74",
}


def region_from_filename(path: Path) -> str | None:
    stable_match = re.fullmatch(
        r"(PMD|LPB|VHPC|PVH|BLA|DHPC|CEA|MSC|SSC|IL|PL)_protein_groups\.csv",
        path.name,
    )
    if stable_match:
        return stable_match.group(1)
    match = re.search(r"_(PMD|LPB|VHPC|PVH|BLA|DHPC|CEA|MSC|SSC|IL|PL)_", path.name)
    return match.group(1) if match else None


def has_region_csvs(path: Path) -> bool:
    if not path.exists():
        return False
    return any(region_from_filename(item) in REGION_CODES for item in path.glob("*.csv"))


def default_data_dir() -> Path:
    candidates = [
        PROJECT_ROOT / "data",
        PROJECT_ROOT / "data" / "data1",
        PROJECT_ROOT.parent / "data" / "data1",
        PROJECT_ROOT.parent / "data",
    ]
    for candidate in candidates:
        if has_region_csvs(candidate):
            return candidate
    return PROJECT_ROOT / "data"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=default_data_dir())
    parser.add_argument("--output-pdf", type=Path, default=PROJECT_ROOT / "figures" / "brain_region_pca.pdf")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "results" / "pca")
    parser.add_argument("--font-family", default="Arial")
    exclusion = parser.add_mutually_exclusive_group()
    exclusion.add_argument("--exclude-sample", default=DEFAULT_EXCLUDED_SAMPLE)
    exclusion.add_argument("--include-all", action="store_true")
    parser.add_argument("--max-protein-missing", type=float, default=0.10)
    return parser.parse_args()


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def clean_sample_name(column: str) -> str:
    name = re.sub(r"^\[\d+\]\s+", "", column)
    return re.sub(r"\.raw\.PG\.Quantity$", "", name)


def read_region_report(path: Path) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    region = region_from_filename(path)
    if region is None:
        raise ValueError(f"Cannot determine region from filename: {path.name}")

    table = pd.read_csv(path, na_values=["NaN"], keep_default_na=True, low_memory=False)
    missing_annotations = [column for column in ANNOTATION_COLUMNS if column not in table.columns]
    if missing_annotations:
        raise ValueError(f"{path.name} lacks columns: {missing_annotations}")

    quantity_columns = [column for column in table.columns if column.endswith(".raw.PG.Quantity")]
    if len(quantity_columns) != 24:
        raise ValueError(f"{path.name}: expected 24 quantity columns, found {len(quantity_columns)}")

    target = table["PG.Organisms"].fillna("").str.contains("Mus musculus", regex=False)
    filtered = table.loc[target, ANNOTATION_COLUMNS + quantity_columns].copy()
    if filtered["PG.ProteinAccessions"].duplicated().any():
        duplicates = filtered.loc[
            filtered["PG.ProteinAccessions"].duplicated(keep=False), "PG.ProteinAccessions"
        ].tolist()
        raise ValueError(f"{path.name}: duplicate protein-group accessions: {duplicates[:5]}")

    quantities = filtered.set_index("PG.ProteinAccessions")[quantity_columns]
    quantities = quantities.apply(pd.to_numeric, errors="coerce")
    quantities.columns = [clean_sample_name(column) for column in quantity_columns]
    sample_meta = pd.DataFrame({"Sample": quantities.columns, "Region": region, "SourceFile": path.name})
    return region, quantities, sample_meta


def collect_region_files(data_dir: Path) -> dict[str, Path]:
    region_to_files: dict[str, list[Path]] = {}
    for path in sorted(data_dir.glob("*.csv")):
        if "_QC_" in path.name:
            continue
        region = region_from_filename(path)
        if region in REGION_CODES:
            region_to_files.setdefault(region, []).append(path)

    missing_regions = sorted(set(REGION_CODES) - set(region_to_files))
    if missing_regions:
        raise FileNotFoundError(f"Missing region reports: {missing_regions}")

    duplicate_regions = {region: files for region, files in region_to_files.items() if len(files) > 1}
    if duplicate_regions:
        details = "; ".join(
            f"{region}: {', '.join(path.name for path in files)}"
            for region, files in sorted(duplicate_regions.items())
        )
        raise ValueError(f"Multiple CSV reports matched the same region: {details}")

    return {region: files[0] for region, files in region_to_files.items()}


def add_confidence_ellipse(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
    confidence_scale: float = 2.4477,
) -> None:
    if len(x) < 3:
        return
    covariance = np.cov(x, y)
    if not np.isfinite(covariance).all():
        return
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0)
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    width, height = 2 * confidence_scale * np.sqrt(eigenvalues)
    ax.add_patch(
        Ellipse(
            (np.mean(x), np.mean(y)),
            width=width,
            height=height,
            angle=angle,
            facecolor=color,
            edgecolor=color,
            linewidth=1.2,
            alpha=0.10,
            zorder=1,
        )
    )


def add_region_cluster_label(
    ax: plt.Axes,
    region: str,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
    confidence_scale: float = 2.4477,
) -> None:
    if len(x) == 0:
        return
    center_x = float(np.mean(x))
    center_y = float(np.mean(y))
    vertical_radius = 0.0
    if len(x) >= 3:
        covariance = np.cov(x, y)
        if np.isfinite(covariance).all():
            vertical_radius = confidence_scale * float(np.sqrt(max(covariance[1, 1], 0)))
    local_y_range = max(float(np.nanmax(y) - np.nanmin(y)), 1.0)
    label_y = center_y + vertical_radius + 0.10 * local_y_range
    ax.text(
        center_x,
        label_y,
        region,
        ha="center",
        va="center",
        fontsize=5.4,
        fontweight="bold",
        color=color,
        bbox={
            "boxstyle": "round,pad=0.22,rounding_size=0.18",
            "facecolor": "white",
            "edgecolor": color,
            "linewidth": 1.0,
            "alpha": 0.88,
        },
        zorder=5,
    )


def build_pca(args: argparse.Namespace) -> dict:
    if not 0 <= args.max_protein_missing < 1:
        raise ValueError("--max-protein-missing must be in [0, 1).")

    region_files = collect_region_files(args.data_dir)
    matrices: list[pd.DataFrame] = []
    metadata: list[pd.DataFrame] = []
    input_rows: dict[str, int] = {}
    retained_mouse_rows: dict[str, int] = {}

    for region in REGION_CODES:
        path = region_files[region]
        input_rows[region] = sum(1 for _ in path.open("r", encoding="utf-8")) - 1
        _, matrix, sample_meta = read_region_report(path)
        retained_mouse_rows[region] = matrix.shape[0]
        matrices.append(matrix)
        metadata.append(sample_meta)

    combined = pd.concat(matrices, axis=1, join="inner")
    sample_meta = pd.concat(metadata, ignore_index=True)
    if combined.columns.duplicated().any():
        duplicates = combined.columns[combined.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate clean sample names: {duplicates}")

    samples_before_exclusion = combined.shape[1]
    excluded_sample: str | None = None
    if not args.include_all:
        if args.exclude_sample in combined.columns:
            excluded_sample = args.exclude_sample
            combined = combined.drop(columns=args.exclude_sample)
            sample_meta = sample_meta.loc[sample_meta["Sample"] != args.exclude_sample].copy()
        elif args.exclude_sample != DEFAULT_EXCLUDED_SAMPLE:
            raise ValueError(f"Requested excluded sample not found: {args.exclude_sample}")
        else:
            warnings.warn(
                f"Default excluded sample was not found and will be ignored: {DEFAULT_EXCLUDED_SAMPLE}",
                RuntimeWarning,
            )

    nonpositive_count = int((combined <= 0).sum().sum())
    combined = combined.mask(combined <= 0)
    log2_matrix = np.log2(combined)

    common_mouse_proteins = log2_matrix.shape[0]
    missing_fraction = log2_matrix.isna().mean(axis=1)
    filtered = log2_matrix.loc[missing_fraction <= args.max_protein_missing].copy()
    proteins_after_missing_filter = filtered.shape[0]

    protein_medians = filtered.median(axis=1, skipna=True)
    if protein_medians.isna().any():
        raise ValueError("At least one retained protein group has no observed values.")
    imputed_values = int(filtered.isna().sum().sum())
    filtered = filtered.T.fillna(protein_medians).T
    filtered = filtered.loc[filtered.var(axis=1, ddof=1) > 0]
    proteins_after_variance_filter = filtered.shape[0]
    if proteins_after_variance_filter < 2:
        raise ValueError("Too few protein groups remain for PCA.")

    sample_by_protein = filtered.T
    n_components = min(10, sample_by_protein.shape[0], sample_by_protein.shape[1])
    model = PCA(n_components=n_components, svd_solver="full")
    coordinates = model.fit_transform(sample_by_protein)

    scores = pd.DataFrame(
        coordinates,
        index=sample_by_protein.index,
        columns=[f"PC{i}" for i in range(1, n_components + 1)],
    ).reset_index(names="Sample")
    scores = scores.merge(sample_meta, on="Sample", how="left", validate="one_to_one")
    if scores["Region"].isna().any():
        raise ValueError("PCA scores could not be matched to all sample metadata.")

    variance = pd.DataFrame(
        {
            "PC": [f"PC{i}" for i in range(1, n_components + 1)],
            "ExplainedVarianceRatio": model.explained_variance_ratio_,
            "ExplainedVariancePercent": 100 * model.explained_variance_ratio_,
        }
    )

    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output_dir / "brain_region_pca_scores.csv", index=False)
    variance.to_csv(args.output_dir / "brain_region_pca_variance.csv", index=False)

    plt.rcParams.update({
        "font.family": args.font_family,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": True,
    })
    markers = dict(zip(REGION_CODES, ["o", "s", "^", "D", "P", "X", "v", "<", ">", "h", "*"]))
    fig, ax = plt.subplots(figsize=(4.8, 4.0), constrained_layout=False)
    fig.subplots_adjust(left=0.12, right=0.78, bottom=0.12, top=0.92)

    for region in REGION_CODES:
        group = scores.loc[scores["Region"] == region]
        color = REGION_PALETTE[region]
        add_confidence_ellipse(ax, group["PC1"].to_numpy(), group["PC2"].to_numpy(), color)
        ax.scatter(
            group["PC1"],
            group["PC2"],
            s=12,
            marker=markers[region],
            color=color,
            edgecolor="white",
            linewidth=0.24,
            alpha=0.88,
            label=f"{region} (n={len(group)})",
            zorder=3,
        )
        add_region_cluster_label(ax, region, group["PC1"].to_numpy(), group["PC2"].to_numpy(), color)

    ax.axhline(0, color="#C9CDD2", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#C9CDD2", linewidth=0.8, zorder=0)
    ax.grid(color="#E7E9EC", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=5.8)
    ax.set_xlabel(f"PC1 ({variance.loc[0, 'ExplainedVariancePercent']:.1f}%)", fontsize=6.8)
    ax.set_ylabel(f"PC2 ({variance.loc[1, 'ExplainedVariancePercent']:.1f}%)", fontsize=6.8)
    ax.set_title("")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(
        title="Brain region",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        borderaxespad=0,
        labelspacing=0.75,
        fontsize=5.4,
        title_fontsize=5.8,
    )
    fig.savefig(args.output_pdf, format="pdf", dpi=300)
    plt.close(fig)

    region_counts = scores.groupby("Region", sort=False).size().reindex(REGION_CODES).to_dict()
    summary = {
        "input_csv_files": len(region_files),
        "samples_before_exclusion": samples_before_exclusion,
        "excluded_sample": excluded_sample,
        "samples_after_exclusion": int(scores.shape[0]),
        "samples_per_region": {key: int(value) for key, value in region_counts.items()},
        "input_protein_groups_per_region": input_rows,
        "mouse_protein_groups_per_region": retained_mouse_rows,
        "shared_mouse_protein_groups": common_mouse_proteins,
        "max_protein_missing_fraction": args.max_protein_missing,
        "protein_groups_after_missing_filter": proteins_after_missing_filter,
        "protein_groups_after_variance_filter": proteins_after_variance_filter,
        "nonpositive_values_treated_as_missing": nonpositive_count,
        "missing_values_median_imputed": imputed_values,
        "remaining_missing_values": int(filtered.isna().sum().sum()),
        "finite_pca_scores": bool(np.isfinite(coordinates).all()),
        "explained_variance_percent": {
            row.PC: float(row.ExplainedVariancePercent)
            for row in variance.itertuples(index=False)
        },
        "output_pdf": portable_path(args.output_pdf),
    }
    with (args.output_dir / "brain_region_pca_qa.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    return summary


def main() -> None:
    args = parse_args()
    summary = build_pca(args)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
