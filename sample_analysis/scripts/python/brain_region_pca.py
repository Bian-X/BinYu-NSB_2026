#!/usr/bin/env python3
"""Region-specific >=20 detection filter, log2, observed sample-median
normalization, protein-wise median filling, accession inner join and PCA.
Exclusion precedes preprocessing. Default: independent 263/264 variants.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from scipy.stats import chi2


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
OUTLIER_SAMPLE = "Astral_MSP2600441_253_G6_1_PVH"
IRT_STANDARD_ACCESSIONS = {"iRT-Kit_WR_fusion"}
MIN_DETECTED_SAMPLES = 20  # >= 20 of 24 biological samples per region
REGION_PALETTE = {'PMD': '#8c510a', 'LPB': '#636363', 'VHPC': '#a6dba0', 'PVH': '#d73027', 'BLA': '#6a3d9a', 'DHPC': '#1b7837', 'CEA': '#c2a5cf', 'MSC': '#fdb863', 'SSC': '#e66101', 'IL': '#1f78b4', 'PL': '#a6cee3'}


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
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "results" / "pca")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "figures")
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--point-size", type=float, default=5.0)
    parser.add_argument("--include-all", action="store_true",
                        help="Produce only the retained-outlier sensitivity variant (default: both variants)")
    parser.add_argument("--imputation", choices=("row-median",), default="row-median")
    parser.add_argument("--exclusions", type=Path, default=PROJECT_ROOT / "sample_exclusions.csv")
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


def is_irt_standard(accessions: pd.Series) -> pd.Series:
    def hit(value: str) -> bool:
        return any(part.strip() in IRT_STANDARD_ACCESSIONS for part in str(value).split(";"))

    return accessions.fillna("").map(hit)


def read_and_preprocess_region(path: Path, excluded_samples=(), imputation="row-median") -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Reference-workflow preprocessing for one region report."""
    region = region_from_filename(path)
    if region is None:
        raise ValueError(f"Cannot determine region from filename: {path.name}")

    table = pd.read_csv(path, na_values=["NaN"], keep_default_na=True, low_memory=False)
    input_label_mapping = [{"source_column": c, "analysis_column": c} for c in table.columns if c.endswith(".raw.PG.Quantity")]
    if table.columns.duplicated().any():
        raise ValueError(f"{path.name}: duplicate input columns")
    missing_annotations = [column for column in ANNOTATION_COLUMNS if column not in table.columns]
    if missing_annotations:
        raise ValueError(f"{path.name} lacks columns: {missing_annotations}")

    quantity_columns = [column for column in table.columns if column.endswith(".raw.PG.Quantity")]
    if len(quantity_columns) != 24:
        raise ValueError(f"{path.name}: expected 24 quantity columns, found {len(quantity_columns)}")

    identities = []
    for column in quantity_columns:
        sample = clean_sample_name(column)
        match = re.fullmatch(r"Astral_MSP2600441_\d+_(G[1-6])_([1-4])_([A-Z]+)(?:_\d{14})?", sample)
        if not match or match.group(3) != region:
            raise ValueError(f"{path.name}: invalid sample identity or region: {sample}")
        identities.append((match.group(1), match.group(2)))
    if len(set(identities)) != 24 or set(identities) != {(f"G{g}",str(r)) for g in range(1,7) for r in range(1,5)}:
        raise ValueError(f"{path.name}: expected unique G1-G6 x replicates 1-4")
    if len({clean_sample_name(c) for c in quantity_columns}) != 24:
        raise ValueError(f"{path.name}: duplicate sample identifiers")

    # Step 1: mouse proteins only, iRT standard removed.
    mouse = table["PG.Organisms"].fillna("").str.contains("Mus musculus", regex=False)
    keep = mouse & ~is_irt_standard(table["PG.ProteinAccessions"])
    filtered = table.loc[keep, ANNOTATION_COLUMNS + quantity_columns].copy()
    n_after_clean = len(filtered)
    if filtered["PG.ProteinAccessions"].duplicated().any():
        duplicates = filtered.loc[
            filtered["PG.ProteinAccessions"].duplicated(keep=False), "PG.ProteinAccessions"
        ].tolist()
        raise ValueError(f"{path.name}: duplicate protein-group accessions: {duplicates[:5]}")

    quantities = filtered.set_index("PG.ProteinAccessions")[quantity_columns]
    quantities = quantities.apply(pd.to_numeric, errors="coerce")
    quantities = quantities.mask((quantities <= 0) | ~np.isfinite(quantities))
    quantities.columns = [clean_sample_name(column) for column in quantity_columns]

    # Exclude before detection filtering, imputation and normalization.
    to_drop = [sample for sample in excluded_samples if sample.endswith("_" + region)]
    missing_excluded = set(to_drop) - set(quantities.columns)
    if missing_excluded:
        raise ValueError(f"Exclusion policy samples not found: {sorted(missing_excluded)}")
    input_samples = list(quantities.columns)
    quantities = quantities.drop(columns=to_drop)
    if len(quantities.columns) < MIN_DETECTED_SAMPLES:
        raise ValueError(f"{region}: insufficient samples for the fixed >=20 filter")
    # The fixed >=20 count is retained in both variants (20/23 for primary PVH).
    detected = quantities.notna().sum(axis=1)
    quantities = quantities.loc[detected >= MIN_DETECTED_SAMPLES]
    n_after_filter = len(quantities)
    if n_after_filter == 0:
        raise ValueError(f"{path.name}: no proteins passed the {MIN_DETECTED_SAMPLES}/24 filter.")

    # Step 3: normalize observed log2 quantities before any filling.
    if imputation != "row-median":
        raise ValueError("Only row-median imputation is supported")
    imputed_cells = int(quantities.isna().sum().sum())
    log2_matrix = np.log2(quantities)
    sample_medians = log2_matrix.median(axis=0)
    if not np.isfinite(sample_medians).all():
        raise ValueError(f"{path.name}: a retained sample has no observations")
    normalized = log2_matrix.subtract(sample_medians, axis=1).add(sample_medians.median())
    # Step 4: one median per protein group on the normalized log2 scale.
    medians = normalized.median(axis=1)
    normalized = normalized.T.fillna(medians).T
    if not np.isfinite(normalized.to_numpy()).all():
        raise ValueError(f"{path.name}: nonfinite values after median imputation")

    sample_meta = pd.DataFrame({"Sample": normalized.columns, "Region": region, "SourceFile": path.name})
    qc = {
        "input_label_mapping": input_label_mapping,
        "region": region,
        "proteins_after_cleaning": int(n_after_clean),
        "proteins_after_detection_filter": int(n_after_filter),
        "imputed_cells": imputed_cells,
        "imputation": imputation,
        "samples_retained": len(quantities.columns),
        "detection_threshold_count": MIN_DETECTED_SAMPLES,
        "sample_manifest": [{"sample_id": c, "included": c not in to_drop} for c in input_samples],
    }
    return normalized, sample_meta, qc


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
    confidence_scale: float | None = None,
) -> None:
    """Approximate 95% bivariate data ellipse, not a CI for the mean."""
    if len(x) < 3:
        return
    if confidence_scale is None:
        confidence_scale = float(np.sqrt(chi2.ppf(.95, df=2)))
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


def square_pca_axes(ax: plt.Axes) -> None:
    """Square plotting area with equal PC1/PC2 distance units."""
    xmin,xmax=ax.get_xlim();ymin,ymax=ax.get_ylim()
    span=max(xmax-xmin,ymax-ymin)
    xmid=(xmin+xmax)/2;ymid=(ymin+ymax)/2
    ax.set_xlim(xmid-span/2,xmid+span/2)
    ax.set_ylim(ymid-span/2,ymid+span/2)
    ax.set_box_aspect(1)
    ax.set_aspect('equal',adjustable='box')


def plot_pca(scores: pd.DataFrame, variance: pd.DataFrame, output_pdf: Path,
             font_family: str, point_size: float, title_suffix: str) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": True,
    })
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=False)
    fig.subplots_adjust(left=0.12, right=0.76, bottom=0.12, top=0.92)

    outlier_scores = scores.loc[scores["Sample"] == OUTLIER_SAMPLE]

    for region in REGION_CODES:
        group = scores.loc[scores["Region"] == region]
        if group.empty:
            continue
        color = REGION_PALETTE[region]
        add_confidence_ellipse(ax, group["PC1"].to_numpy(), group["PC2"].to_numpy(), color)
        ax.scatter(
            group["PC1"],
            group["PC2"],
            s=point_size,
            marker="o",
            color=color,
            edgecolor="white",
            linewidth=0.2,
            alpha=0.88,
            label=f"{region} (n={len(group)})",
            zorder=3,
        )

    # Mark the low-coverage outlier sample when present.
    if not outlier_scores.empty:
        ax.scatter(
            outlier_scores["PC1"],
            outlier_scores["PC2"],
            s=point_size * 4,
            marker="o",
            facecolor="none",
            edgecolor="black",
            linewidth=1.0,
            label="PVH outlier (G6_1)",
            zorder=4,
        )

    ax.axhline(0, color="#C9CDD2", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#C9CDD2", linewidth=0.8, zorder=0)
    ax.grid(color="#E7E9EC", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=5.8)
    ax.set_xlabel(f"PC1 ({variance.loc[0, 'ExplainedVariancePercent']:.1f}%)", fontsize=6.8)
    ax.set_ylabel(f"PC2 ({variance.loc[1, 'ExplainedVariancePercent']:.1f}%)", fontsize=6.8)
    ax.set_title(title_suffix, fontsize=7.2)
    square_pca_axes(ax)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(
        title="Brain region",
        loc="upper left",
        bbox_to_anchor=(1.04, 1),
        ncol=1,
        frameon=False,
        borderaxespad=0,
        labelspacing=0.75,
        fontsize=5.4,
        title_fontsize=5.8,
    )
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, format="pdf", dpi=300)
    fig.savefig(output_pdf.with_suffix('.png'), dpi=400)
    plt.close(fig)


def run_pca(matrix: pd.DataFrame, sample_meta: pd.DataFrame, label: str,
            args: argparse.Namespace) -> dict:
    sample_by_protein = matrix.T
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

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output_dir / f"brain_region_pca_{label}_scores.csv", index=False)
    variance.to_csv(args.output_dir / f"brain_region_pca_{label}_variance.csv", index=False)

    title_suffix = (
        f"{len(scores)} samples (PVH outlier retained)"
        if label == "include_all"
        else f"{len(scores)} samples (PVH outlier excluded before preprocessing)"
    )
    output_pdf = args.figure_dir / f"brain_region_pca_{label}.pdf"
    plot_pca(scores, variance, output_pdf, args.font_family, args.point_size, title_suffix)

    return {
        "samples": int(scores.shape[0]),
        "explained_variance_percent": {
            row.PC: float(row.ExplainedVariancePercent)
            for row in variance.itertuples(index=False)
        },
        "output_pdf": portable_path(output_pdf),
    }


def main() -> None:
    args = parse_args()

    region_files = collect_region_files(args.data_dir)
    policy = pd.read_csv(args.exclusions, keep_default_na=False)
    if not {"sample_id", "reason"}.issubset(policy.columns) or policy["sample_id"].duplicated().any() or policy["sample_id"].eq("").any():
        raise ValueError("Invalid sample exclusion policy")
    exclusions = policy["sample_id"].tolist()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "preprocessing": "per region; exclusion before fixed >=20 detection count; " + args.imputation +
                         " imputation on normalized log2 scale; accession inner join",
        "outlier_sample": OUTLIER_SAMPLE,
        "variants": {},
    }
    labels = ("include_all",) if args.include_all else ("exclude_outlier", "include_all")
    for label in labels:
        matrices, metadata, qc_records = [], [], []
        for region in REGION_CODES:
            normalized, sample_meta, qc = read_and_preprocess_region(
                region_files[region], exclusions if label == "exclude_outlier" else (), args.imputation)
            matrices.append(normalized)
            metadata.append(sample_meta)
            qc_records.append(qc)
        combined = pd.concat(matrices, axis=1, join="inner")
        sample_meta = pd.concat(metadata, ignore_index=True)
        if combined.columns.duplicated().any():
            raise ValueError("Duplicate clean sample names")
        variant = run_pca(combined, sample_meta, label, args)
        variant["shared_protein_groups_after_inner_join"] = int(combined.shape[0])
        variant["per_region_qc"] = qc_records
        summary["variants"][label] = variant
        pd.DataFrame([entry for qc in qc_records for entry in qc["sample_manifest"]]).to_csv(
            args.output_dir / f"brain_region_pca_{label}_sample_manifest.csv", index=False)
    with (args.output_dir / "brain_region_pca_qa.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
