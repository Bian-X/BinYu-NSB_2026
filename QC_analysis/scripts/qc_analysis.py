#!/usr/bin/env python
"""Compute reproducible QC summaries from DIA-LFQ protein-group reports.

The PDF entry point imports the data-loading and calculation functions in this
module. The legacy SVG command-line entry point is retained only to reproduce
the originally archived supplementary SVG files; use ``run_qc_pdf.py`` for the
GitHub release workflow.
"""

from __future__ import annotations

import argparse
import html
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ANNOTATION_COLUMNS = [
    "PG.ProteinAccessions",
    "PG.Genes",
    "PG.Organisms",
    "PG.ProteinDescriptions",
    "PG.NrOfStrippedSequencesIdentified (Experiment-wide)",
]

REGION_ORDER = [
    "QC",
    "PL",
    "IL",
    "SSC",
    "MSC",
    "BLA",
    "CEA",
    "DHPC",
    "PVH",
    "VHPC",
    "PMD",
    "LPB",
]

GROUP_ORDER = ["QC", "G1", "G2", "G3", "G4", "G5", "G6"]

GROUP_DEFINITIONS = {
    "QC": ("pooled_quality_control", "technical_qc"),
    "G1": ("untreated_home_cage_control", "baseline_control"),
    "G2": ("temporally_unpaired_tone_shock_control", "non_associative_control"),
    "G3": ("single_6s_shock_only_context_learning", "learning_control"),
    "G4": ("paired_tone_shock_learning", "learning"),
    "G5": ("contextual_memory_retrieval", "retrieval"),
    "G6": ("tone_cued_memory_retrieval", "retrieval"),
}

REGION_COLORS = {
    "QC": "#2f3437",
    "BLA": "#6a3d9a",
    "CEA": "#c2a5cf",
    "IL": "#1f78b4",
    "PL": "#a6cee3",
    "DHPC": "#1b7837",
    "VHPC": "#a6dba0",
    "SSC": "#e66101",
    "MSC": "#fdb863",
    "PVH": "#d73027",
    "PMD": "#8c510a",
    "LPB": "#636363",
    "UNKNOWN": "#8c8c8c",
}

# Figure 2F groups anatomically related regions into shared color families.
HEATMAP_REGION_ORDER = ["QC", "BLA", "CEA", "IL", "PL", "DHPC", "VHPC", "SSC", "MSC", "PVH", "PMD", "LPB"]
# Heatmap annotations use the global brain-region palette used by all QC figures.
HEATMAP_REGION_COLORS = REGION_COLORS


@dataclass
class ReportData:
    file_path: Path
    region: str
    is_qc_report: bool
    annotation: pd.DataFrame
    quantity: pd.DataFrame
    sample_metadata: pd.DataFrame


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def clean_quantity_column(col: str) -> str:
    cleaned = re.sub(r"^\[\d+\]\s*", "", col)
    cleaned = cleaned.replace(".PG.Quantity", "")
    return cleaned


def parse_report_region(path: Path) -> str:
    match = re.search(r"Astral_MSP2600441_([^_]+)_DIA", path.name)
    if match:
        return match.group(1)
    return "UNKNOWN"


def parse_sample(raw_name: str, report_file: str, report_region: str) -> dict[str, object]:
    qc_match = re.match(r"^Astral_MSP2600441_QC_(\d+)\.raw$", raw_name)
    if qc_match:
        qc_index = int(qc_match.group(1))
        condition, level = GROUP_DEFINITIONS["QC"]
        return {
            "sample_id": raw_name.replace(".raw", ""),
            "raw_file": raw_name,
            "report_file": report_file,
            "sample_type": "QC",
            "region": "QC",
            "group": "QC",
            "condition": condition,
            "level": level,
            "replicate": qc_index,
            "filename_index": np.nan,
            "qc_index": qc_index,
            "timestamp_suffix": "",
            "report_region": report_region,
        }

    pattern = (
        r"^Astral_MSP2600441_(\d+)_(G[1-6])_(\d+)_"
        r"([A-Za-z0-9]+)(?:_(\d{14}))?\.raw$"
    )
    match = re.match(pattern, raw_name)
    if match:
        filename_index, group, replicate, region, timestamp = match.groups()
        condition, level = GROUP_DEFINITIONS[group]
        return {
            "sample_id": raw_name.replace(".raw", ""),
            "raw_file": raw_name,
            "report_file": report_file,
            "sample_type": "BIO",
            "region": region,
            "group": group,
            "condition": condition,
            "level": level,
            "replicate": int(replicate),
            "filename_index": int(filename_index),
            "qc_index": np.nan,
            "timestamp_suffix": timestamp or "",
            "report_region": report_region,
        }

    return {
        "sample_id": raw_name.replace(".raw", ""),
        "raw_file": raw_name,
        "report_file": report_file,
        "sample_type": "UNKNOWN",
        "region": report_region,
        "group": "UNKNOWN",
        "condition": "UNKNOWN",
        "level": "",
        "replicate": np.nan,
        "filename_index": np.nan,
        "qc_index": np.nan,
        "timestamp_suffix": "",
        "report_region": report_region,
    }


def region_rank(region: str) -> int:
    return REGION_ORDER.index(region) if region in REGION_ORDER else len(REGION_ORDER)


def group_rank(group: str) -> int:
    return GROUP_ORDER.index(group) if group in GROUP_ORDER else len(GROUP_ORDER)


def add_analysis_decisions(table: pd.DataFrame) -> pd.DataFrame:
    """Add rerun and downstream-analysis decisions to a sample-level table."""
    result = table.copy()
    result["is_rerun"] = np.where(
        (result["sample_type"] == "BIO") & result["timestamp_suffix"].fillna("").ne(""),
        "yes",
        "no",
    )
    result["analysis_inclusion"] = "yes"
    result["qc_flag"] = np.where(result["sample_type"] == "QC", "technical_qc", "retain")
    result["exclusion_reason"] = ""

    pvh_sample = "Astral_MSP2600441_253_G6_1_PVH"
    cea_sample = "Astral_MSP2600441_112_G3_4_CEA"

    pvh = result["sample_id"].eq(pvh_sample)
    result.loc[pvh, "analysis_inclusion"] = "no"
    result.loc[pvh, "qc_flag"] = "exclude_extreme_low_coverage"
    result.loc[pvh, "exclusion_reason"] = (
        "33.21% sample-level missingness; 6,143 detected protein groups"
    )

    cea = result["sample_id"].eq(cea_sample)
    result.loc[cea, "qc_flag"] = "retain_elevated_missingness"
    result.loc[cea, "exclusion_reason"] = (
        "10.87% sample-level missingness; retained for analysis"
    )
    return result


def load_report(path: Path) -> ReportData:
    region = parse_report_region(path)
    is_qc_report = region == "QC"
    df = pd.read_csv(path, na_values=["NaN", "nan", ""], keep_default_na=True)
    quantity_cols = [c for c in df.columns if c.endswith(".PG.Quantity")]
    if "PG.ProteinAccessions" not in df.columns:
        raise ValueError(f"{path.name} does not contain the required PG.ProteinAccessions column.")
    if not quantity_cols:
        raise ValueError(f"{path.name} does not contain any *.PG.Quantity columns.")
    annotation_cols = [c for c in ANNOTATION_COLUMNS if c in df.columns]
    annotation = df[annotation_cols].copy()

    sample_rows = []
    rename_map = {}
    for col in quantity_cols:
        raw_name = clean_quantity_column(col)
        parsed = parse_sample(raw_name, path.name, region)
        sample_rows.append(parsed)
        rename_map[col] = parsed["sample_id"]

    quantity = df[quantity_cols].rename(columns=rename_map)
    quantity.index = annotation["PG.ProteinAccessions"].astype(str)
    quantity = quantity.apply(pd.to_numeric, errors="coerce")
    if quantity.index.duplicated().any():
        count = int(quantity.index.duplicated().sum())
        raise ValueError(f"{path.name} contains {count} duplicated protein-group accession(s).")
    if quantity.columns.duplicated().any():
        count = int(quantity.columns.duplicated().sum())
        raise ValueError(f"{path.name} contains {count} duplicated sample identifier(s).")
    sample_metadata = pd.DataFrame(sample_rows)
    return ReportData(path, region, is_qc_report, annotation, quantity, sample_metadata)


def load_all_reports(data_dir: Path) -> list[ReportData]:
    files = sorted(data_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found under {data_dir}")
    return [load_report(path) for path in files]


def build_combined_matrix(reports: list[ReportData]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_sample_meta = pd.concat([r.sample_metadata for r in reports], ignore_index=True)
    if all_sample_meta["sample_id"].duplicated().any():
        duplicates = all_sample_meta.loc[all_sample_meta["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Duplicate sample identifier(s) across reports: {duplicates[:5]}")

    annotation_frames = []
    matrix_frames = []
    for report in reports:
        ann = report.annotation.copy()
        ann["PG.ProteinAccessions"] = ann["PG.ProteinAccessions"].astype(str)
        annotation_frames.append(ann)
        matrix_frames.append(report.quantity)

    annotation = pd.concat(annotation_frames, ignore_index=True)
    annotation = annotation.drop_duplicates("PG.ProteinAccessions", keep="first")
    annotation = annotation.set_index("PG.ProteinAccessions", drop=False)

    matrix = pd.concat(matrix_frames, axis=1, join="outer")
    matrix = matrix.groupby(matrix.index).mean(numeric_only=True)

    sort_cols = ["region_rank", "group_rank", "replicate", "filename_index", "qc_index", "sample_id"]
    all_sample_meta["region_rank"] = all_sample_meta["region"].map(region_rank)
    all_sample_meta["group_rank"] = all_sample_meta["group"].map(group_rank)
    all_sample_meta = all_sample_meta.sort_values(sort_cols, na_position="last").reset_index(drop=True)
    matrix = matrix[[c for c in all_sample_meta["sample_id"] if c in matrix.columns]]
    return matrix, all_sample_meta, annotation


def sample_qc_summary_from_reports(reports: list[ReportData]) -> pd.DataFrame:
    """Compute sample QC metrics inside each report's own protein universe.

    This avoids inflating missingness by using the cross-file union of proteins.
    """

    rows = []
    for report in reports:
        log2 = np.log2(report.quantity)
        for _, meta_row in report.sample_metadata.iterrows():
            sample_id = meta_row["sample_id"]
            values = report.quantity[sample_id]
            log_values = log2[sample_id]
            rows.append(
                {
                    **meta_row.to_dict(),
                    "report_protein_groups": report.quantity.shape[0],
                    "detected_protein_groups": int(values.notna().sum()),
                    "missing_percent": float(values.isna().mean() * 100),
                    "total_intensity": float(values.sum(skipna=True)),
                    "median_intensity": float(values.median(skipna=True)),
                    "median_log2_intensity": float(log_values.median(skipna=True)),
                    "q25_log2_intensity": float(log_values.quantile(0.25)),
                    "q75_log2_intensity": float(log_values.quantile(0.75)),
                }
            )
    qc = pd.DataFrame(rows)
    qc["region_rank"] = qc["region"].map(region_rank)
    qc["group_rank"] = qc["group"].map(group_rank)
    return qc.sort_values(["region_rank", "group_rank", "replicate", "filename_index", "qc_index"]).reset_index(drop=True)


def median_normalize_log2(log2_matrix: pd.DataFrame) -> pd.DataFrame:
    sample_medians = log2_matrix.median(axis=0, skipna=True)
    global_median = float(sample_medians.median(skipna=True))
    return log2_matrix.subtract(sample_medians, axis=1).add(global_median)


def prepare_complete_matrix(matrix: pd.DataFrame, min_detection_fraction: float = 0.8) -> pd.DataFrame:
    if not 0 < min_detection_fraction <= 1:
        raise ValueError("min_detection_fraction must be greater than 0 and no more than 1.")
    log2 = np.log2(matrix)
    norm = median_normalize_log2(log2)
    keep = norm.notna().mean(axis=1) >= min_detection_fraction
    filtered = norm.loc[keep].copy()
    row_medians = filtered.median(axis=1, skipna=True)
    filtered = filtered.T.fillna(row_medians).T
    return filtered


def compute_pca(matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    x = matrix.T.to_numpy(dtype=float)
    x = x - np.nanmean(x, axis=0, keepdims=True)
    feature_sd = np.nanstd(x, axis=0, ddof=1)
    feature_sd[feature_sd == 0] = 1.0
    x = x / feature_sd
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    _u, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = x @ vt[:2].T
    explained = (s**2) / np.sum(s**2)
    result = pd.DataFrame(scores, index=matrix.columns, columns=["PC1", "PC2"])
    return result, explained[:2]


def pearson_corr_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix.corr(method="pearson")


def qc_stability_metrics(reports: list[ReportData]) -> pd.DataFrame:
    qc_reports = [report for report in reports if report.is_qc_report]
    if not qc_reports:
        raise ValueError("No QC report found.")
    report = qc_reports[0]
    sample_meta = report.sample_metadata.copy()
    sample_meta["region_rank"] = sample_meta["region"].map(region_rank)
    sample_meta["group_rank"] = sample_meta["group"].map(group_rank)
    qc_samples = sample_meta["sample_id"].tolist()
    qc_matrix = report.quantity[qc_samples]
    log2_qc = np.log2(qc_matrix)
    median_profile = log2_qc.median(axis=1, skipna=True)
    rows = []
    for sample_id in qc_samples:
        values = log2_qc[sample_id]
        valid = values.notna() & median_profile.notna()
        if valid.sum() >= 3:
            corr = float(np.corrcoef(values[valid], median_profile[valid])[0, 1])
        else:
            corr = np.nan
        meta = sample_meta.loc[sample_meta["sample_id"] == sample_id].iloc[0].to_dict()
        rows.append(
            {
                **meta,
                "detected_protein_groups": int(qc_matrix[sample_id].notna().sum()),
                "missing_percent": float(qc_matrix[sample_id].isna().mean() * 100),
                "median_log2_intensity": float(values.median(skipna=True)),
                "correlation_to_qc_median": corr,
            }
        )
    return pd.DataFrame(rows).sort_values("qc_index")


def qc_cv_table(matrix: pd.DataFrame, sample_meta: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    qc_samples = sample_meta.loc[sample_meta["sample_type"] == "QC", "sample_id"].tolist()
    qc_matrix = matrix[qc_samples]
    mean_intensity = qc_matrix.mean(axis=1, skipna=True)
    sd_intensity = qc_matrix.std(axis=1, skipna=True)
    detected = qc_matrix.notna().sum(axis=1)
    cv_percent = (sd_intensity / mean_intensity) * 100
    log2_qc = np.log2(qc_matrix)
    table = pd.DataFrame(
        {
            "PG.ProteinAccessions": qc_matrix.index,
            "qc_detected_n": detected.values,
            "qc_detected_fraction": (detected / len(qc_samples)).values,
            "qc_mean_intensity": mean_intensity.values,
            "qc_mean_log2_intensity": log2_qc.mean(axis=1, skipna=True).values,
            "qc_sd_log2_intensity": log2_qc.std(axis=1, skipna=True).values,
            "qc_cv_percent": cv_percent.values,
        }
    )
    ann_cols = [c for c in ANNOTATION_COLUMNS if c in annotation.columns and c != "PG.ProteinAccessions"]
    if ann_cols:
        ann_for_merge = annotation.reset_index(drop=True)[["PG.ProteinAccessions"] + ann_cols]
        table = table.merge(
            ann_for_merge,
            on="PG.ProteinAccessions",
            how="left",
        )
    return table.sort_values(["qc_detected_fraction", "qc_cv_percent"], ascending=[False, True])


def protein_detection_by_file(reports: list[ReportData]) -> pd.DataFrame:
    rows = []
    for report in reports:
        detected_per_protein = report.quantity.notna().sum(axis=1)
        rows.append(
            {
                "report_file": report.file_path.name,
                "region": report.region,
                "sample_columns": report.quantity.shape[1],
                "protein_groups": report.quantity.shape[0],
                "missing_percent": float(report.quantity.isna().mean().mean() * 100),
                "fully_detected_protein_groups": int((detected_per_protein == report.quantity.shape[1]).sum()),
                "detected_in_at_least_half_samples": int(
                    (detected_per_protein >= math.ceil(report.quantity.shape[1] * 0.5)).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def scale(value: float, old_min: float, old_max: float, new_min: float, new_max: float) -> float:
    if not np.isfinite(value):
        return (new_min + new_max) / 2
    if old_max == old_min:
        return (new_min + new_max) / 2
    return new_min + (value - old_min) * (new_max - new_min) / (old_max - old_min)


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#222}",
        ".title{font-size:24px;font-weight:700}",
        ".subtitle{font-size:13px;fill:#555}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#ddd;stroke-width:1}",
        ".tick{font-size:10px;fill:#555}",
        ".label{font-size:12px;fill:#333}",
        ".legend{font-size:11px;fill:#333}",
        "</style>",
    ]


def save_svg(parts: list[str], path: Path) -> None:
    path.write_text("\n".join(parts + ["</svg>\n"]), encoding="utf-8")


def add_text(parts: list[str], x: float, y: float, text: object, cls: str = "label", anchor: str = "start") -> None:
    parts.append(f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{esc(text)}</text>')


def add_axes(parts: list[str], x0: float, y0: float, w: float, h: float) -> None:
    parts.append(f'<line x1="{x0}" y1="{y0+h}" x2="{x0+w}" y2="{y0+h}" class="axis"/>')
    parts.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+h}" class="axis"/>')


def add_y_ticks(
    parts: list[str],
    x0: float,
    y0: float,
    h: float,
    vmin: float,
    vmax: float,
    ticks: int = 5,
    fmt: str = "{:.0f}",
) -> None:
    for i in range(ticks):
        value = vmin + (vmax - vmin) * i / (ticks - 1)
        y = scale(value, vmin, vmax, y0 + h, y0)
        parts.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0-5}" y2="{y:.1f}" class="axis"/>')
        parts.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+960}" y2="{y:.1f}" class="grid" opacity="0.45"/>')
        add_text(parts, x0 - 8, y + 4, fmt.format(value), "tick", "end")


def add_region_legend(parts: list[str], x: float, y: float, regions: list[str]) -> None:
    for i, region in enumerate(regions):
        yy = y + i * 18
        color = REGION_COLORS.get(region, REGION_COLORS["UNKNOWN"])
        parts.append(f'<rect x="{x}" y="{yy-10}" width="11" height="11" fill="{color}"/>')
        add_text(parts, x + 16, yy, region, "legend")


def figure_1_sample_overview(sample_qc: pd.DataFrame, out_path: Path) -> None:
    df = sample_qc.sort_values(["region_rank", "group_rank", "replicate", "filename_index", "qc_index"]).reset_index(drop=True)
    width, height = 1800, 820
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 2D. Sample detection depth and missingness overview", "title")
    add_text(parts, 40, 60, "Bars show detected protein groups; points show sample missing percentage.", "subtitle")

    x0, y0, w, h = 80, 95, 1600, 285
    add_axes(parts, x0, y0, w, h)
    det_min, det_max = 0, float(df["detected_protein_groups"].max() * 1.05)
    add_y_ticks(parts, x0, y0, h, det_min, det_max, ticks=5)
    n = len(df)
    gap = 1.0
    bar_w = max(1.0, (w - gap * (n - 1)) / n)
    x1, y1, h2 = x0, 470, 250
    region_boundaries = [
        index for index in range(1, n) if df.loc[index, "region"] != df.loc[index - 1, "region"]
    ]
    for index in region_boundaries:
        x = x0 + index * (bar_w + gap) - gap / 2
        parts.append(
            f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y1+h2}" '
            'stroke="#9e9e9e" stroke-dasharray="4,4" stroke-width="1" opacity="0.75"/>'
        )
    for i, row in df.iterrows():
        x = x0 + i * (bar_w + gap)
        y = scale(row["detected_protein_groups"], det_min, det_max, y0 + h, y0)
        color = REGION_COLORS.get(row["region"], REGION_COLORS["UNKNOWN"])
        opacity = "0.95" if row["sample_type"] == "QC" else "0.70"
        parts.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{y0+h-y:.2f}" '
            f'fill="{color}" opacity="{opacity}"><title>{esc(row["sample_id"])}: '
            f'{int(row["detected_protein_groups"])} detected</title></rect>'
        )
    add_text(parts, 20, y0 + h / 2, "Detected protein groups", "label")

    add_axes(parts, x1, y1, w, h2)
    miss_min, miss_max = 0, max(10.0, float(df["missing_percent"].max() * 1.1))
    add_y_ticks(parts, x1, y1, h2, miss_min, miss_max, ticks=6, fmt="{:.0f}%")
    for i, row in df.iterrows():
        x = x1 + i * (bar_w + gap) + bar_w / 2
        y = scale(row["missing_percent"], miss_min, miss_max, y1 + h2, y1)
        color = REGION_COLORS.get(row["region"], REGION_COLORS["UNKNOWN"])
        radius = 3.5 if row["sample_type"] == "QC" else 2.2
        parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{color}" opacity="0.85">'
            f'<title>{esc(row["sample_id"])}: {row["missing_percent"]:.2f}% missing</title></circle>'
        )
    threshold_y = scale(10, miss_min, miss_max, y1 + h2, y1)
    parts.append(
        f'<line x1="{x1}" y1="{threshold_y:.1f}" x2="{x1+w}" y2="{threshold_y:.1f}" '
        'stroke="#b2182b" stroke-dasharray="6,5" stroke-width="1.5"/>'
    )
    add_text(parts, x1 + w + 8, threshold_y + 4, "10% missing", "tick")
    add_text(parts, 20, y1 + h2 / 2, "Missing rate", "label")

    current = None
    start_i = 0
    y_label = 750
    for i, region in enumerate(df["region"].tolist() + ["END"]):
        if current is None:
            current = region
            start_i = i
        elif region != current:
            end_i = i - 1
            cx = x0 + ((start_i + end_i) / 2) * (bar_w + gap) + bar_w / 2
            parts.append(f'<line x1="{cx:.1f}" y1="{y1+h2+6}" x2="{cx:.1f}" y2="{y1+h2+16}" class="axis"/>')
            add_text(parts, cx, y_label, current, "tick", "middle")
            current = region
            start_i = i

    legend_regions = [r for r in REGION_ORDER if r in set(df["region"])]
    add_region_legend(parts, 1695, 110, legend_regions)
    save_svg(parts, out_path)


def figure_2_qc_injection(qc_metrics: pd.DataFrame, out_path: Path) -> None:
    df = qc_metrics.sort_values("qc_index").reset_index(drop=True)
    width, height = 1500, 920
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 2. QC injection-order stability", "title")
    add_text(parts, 40, 60, "The 25 QC samples are ordered by QC injection index.", "subtitle")

    panels = [
        ("detected_protein_groups", "Detected protein groups", "{:.0f}"),
        ("missing_percent", "Missing rate (%)", "{:.1f}"),
        ("median_log2_intensity", "Median log2 intensity", "{:.1f}"),
        ("correlation_to_qc_median", "Correlation to QC median", "{:.3f}"),
    ]
    panel_w, panel_h = 610, 290
    positions = [(95, 105), (805, 105), (95, 510), (805, 510)]
    for (col, title, fmt), (x0, y0) in zip(panels, positions):
        values = df[col].astype(float)
        vmin = float(values.min())
        vmax = float(values.max())
        pad = (vmax - vmin) * 0.15 if vmax > vmin else max(1.0, abs(vmax) * 0.05)
        if col == "missing_percent":
            vmin = 0
        else:
            vmin -= pad
        vmax += pad
        add_text(parts, x0, y0 - 18, title, "label")
        add_axes(parts, x0, y0, panel_w, panel_h)
        add_y_ticks(parts, x0, y0, panel_h, vmin, vmax, ticks=5, fmt=fmt)
        points = []
        for _, row in df.iterrows():
            x = scale(float(row["qc_index"]), 1, 25, x0, x0 + panel_w)
            y = scale(float(row[col]), vmin, vmax, y0 + panel_h, y0)
            points.append((x, y))
        path_d = " ".join(
            [f'{"M" if i == 0 else "L"} {x:.1f} {y:.1f}' for i, (x, y) in enumerate(points)]
        )
        parts.append(f'<path d="{path_d}" fill="none" stroke="#2f3437" stroke-width="2"/>')
        median_y = scale(float(values.median()), vmin, vmax, y0 + panel_h, y0)
        parts.append(
            f'<line x1="{x0}" y1="{median_y:.1f}" x2="{x0+panel_w}" y2="{median_y:.1f}" '
            'stroke="#bdbdbd" stroke-dasharray="5,5"/>'
        )
        for _, row in df.iterrows():
            x = scale(float(row["qc_index"]), 1, 25, x0, x0 + panel_w)
            y = scale(float(row[col]), vmin, vmax, y0 + panel_h, y0)
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="#2f3437">'
                f'<title>QC_{int(row["qc_index"])}: {row[col]:.4g}</title></circle>'
            )
        for tick in [1, 5, 10, 15, 20, 25]:
            x = scale(tick, 1, 25, x0, x0 + panel_w)
            parts.append(f'<line x1="{x:.1f}" y1="{y0+panel_h}" x2="{x:.1f}" y2="{y0+panel_h+5}" class="axis"/>')
            add_text(parts, x, y0 + panel_h + 20, tick, "tick", "middle")
        add_text(parts, x0 + panel_w / 2, y0 + panel_h + 42, "QC injection order", "label", "middle")
    save_svg(parts, out_path)


def figure_3_intensity_distribution(sample_qc: pd.DataFrame, log2_matrix: pd.DataFrame, out_path: Path) -> None:
    norm = median_normalize_log2(log2_matrix)
    before = sample_qc.set_index("sample_id")[["region", "sample_type", "region_rank", "group_rank", "replicate"]].copy()
    before["median"] = log2_matrix.median(axis=0)
    before["q25"] = log2_matrix.quantile(0.25, axis=0)
    before["q75"] = log2_matrix.quantile(0.75, axis=0)
    before["panel"] = "Before median normalization"
    after = sample_qc.set_index("sample_id")[["region", "sample_type", "region_rank", "group_rank", "replicate"]].copy()
    after["median"] = norm.median(axis=0)
    after["q25"] = norm.quantile(0.25, axis=0)
    after["q75"] = norm.quantile(0.75, axis=0)
    after["panel"] = "After median normalization"

    width, height = 1800, 760
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 3. Log2 intensity distribution before and after normalization", "title")
    add_text(parts, 40, 60, "Each vertical segment is one sample IQR; the dot is the sample median.", "subtitle")

    combined = pd.concat([before, after])
    vmin = float(combined["q25"].min() - 0.5)
    vmax = float(combined["q75"].max() + 0.5)
    ordered_samples = sample_qc.sort_values(["region_rank", "group_rank", "replicate", "filename_index", "qc_index"])[
        "sample_id"
    ].tolist()
    panel_specs = [(before, 90, 115, "Before median normalization"), (after, 90, 435, "After median normalization")]
    w, h = 1600, 230
    for data, x0, y0, title in panel_specs:
        add_text(parts, x0, y0 - 18, title, "label")
        add_axes(parts, x0, y0, w, h)
        add_y_ticks(parts, x0, y0, h, vmin, vmax, ticks=5, fmt="{:.1f}")
        n = len(ordered_samples)
        for i, sample in enumerate(ordered_samples):
            row = data.loc[sample]
            x = scale(i, 0, n - 1, x0 + 2, x0 + w - 2)
            y25 = scale(row["q25"], vmin, vmax, y0 + h, y0)
            y75 = scale(row["q75"], vmin, vmax, y0 + h, y0)
            ymed = scale(row["median"], vmin, vmax, y0 + h, y0)
            color = REGION_COLORS.get(row["region"], REGION_COLORS["UNKNOWN"])
            parts.append(
                f'<line x1="{x:.1f}" y1="{y25:.1f}" x2="{x:.1f}" y2="{y75:.1f}" '
                f'stroke="{color}" stroke-width="1.1" opacity="0.55"/>'
            )
            radius = 2.8 if row["sample_type"] == "QC" else 1.8
            parts.append(
                f'<circle cx="{x:.1f}" cy="{ymed:.1f}" r="{radius}" fill="{color}" opacity="0.85">'
                f'<title>{esc(sample)} median={row["median"]:.2f}</title></circle>'
            )
    add_region_legend(parts, 1700, 120, [r for r in REGION_ORDER if r in set(sample_qc["region"])])
    save_svg(parts, out_path)


def figure_4_pca(pca_df: pd.DataFrame, explained: np.ndarray, sample_qc: pd.DataFrame, out_path: Path) -> None:
    df = sample_qc.merge(pca_df, left_on="sample_id", right_index=True, how="inner")
    width, height = 1100, 850
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 4. PCA of biological samples and QC samples", "title")
    add_text(parts, 40, 60, "PCA uses median-normalized log2 protein-group quantities with row-median imputation.", "subtitle")
    x0, y0, w, h = 105, 105, 760, 620
    xmin, xmax = float(df["PC1"].min()), float(df["PC1"].max())
    ymin, ymax = float(df["PC2"].min()), float(df["PC2"].max())
    xpad = (xmax - xmin) * 0.08 if xmax > xmin else 1
    ypad = (ymax - ymin) * 0.08 if ymax > ymin else 1
    xmin, xmax = xmin - xpad, xmax + xpad
    ymin, ymax = ymin - ypad, ymax + ypad
    add_axes(parts, x0, y0, w, h)
    add_y_ticks(parts, x0, y0, h, ymin, ymax, ticks=6, fmt="{:.1f}")
    for i in range(6):
        value = xmin + (xmax - xmin) * i / 5
        x = scale(value, xmin, xmax, x0, x0 + w)
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h}" x2="{x:.1f}" y2="{y0+h+5}" class="axis"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0+h}" class="grid" opacity="0.45"/>')
        add_text(parts, x, y0 + h + 20, f"{value:.1f}", "tick", "middle")

    for _, row in df.iterrows():
        x = scale(row["PC1"], xmin, xmax, x0, x0 + w)
        y = scale(row["PC2"], ymin, ymax, y0 + h, y0)
        color = REGION_COLORS.get(row["region"], REGION_COLORS["UNKNOWN"])
        if row["sample_type"] == "QC":
            parts.append(
                f'<rect x="{x-4:.1f}" y="{y-4:.1f}" width="8" height="8" fill="{color}" opacity="0.9">'
                f'<title>{esc(row["sample_id"])}</title></rect>'
            )
        else:
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{color}" opacity="0.75">'
                f'<title>{esc(row["sample_id"])} {esc(row["group"])}</title></circle>'
            )
    add_text(parts, x0 + w / 2, y0 + h + 50, f"PC1 ({explained[0] * 100:.1f}%)", "label", "middle")
    add_text(parts, 28, y0 + h / 2, f"PC2 ({explained[1] * 100:.1f}%)", "label")
    add_region_legend(parts, 900, 120, [r for r in REGION_ORDER if r in set(df["region"])])
    save_svg(parts, out_path)


def color_gradient(value: float, vmin: float = 0.85, vmax: float = 1.0) -> str:
    value = max(vmin, min(vmax, value))
    t = (value - vmin) / (vmax - vmin) if vmax > vmin else 1.0
    # Higher correlation is darker; lower correlation is lighter.
    r = int(222 * (1 - t) + 8 * t)
    g = int(235 * (1 - t) + 69 * t)
    b = int(247 * (1 - t) + 148 * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def figure_5_correlation_heatmap(corr: pd.DataFrame, sample_qc: pd.DataFrame, out_path: Path) -> None:
    heatmap_region_rank = {region: index for index, region in enumerate(HEATMAP_REGION_ORDER)}
    ordered = sample_qc.assign(
        heatmap_region_rank=sample_qc["region"].map(heatmap_region_rank).fillna(len(HEATMAP_REGION_ORDER))
    ).sort_values(["heatmap_region_rank", "group_rank", "replicate", "filename_index", "qc_index"])[
        "sample_id"
    ].tolist()
    ordered = [s for s in ordered if s in corr.index]
    corr = corr.loc[ordered, ordered]
    n = len(ordered)
    width, height = 1500, 1340
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 2F. Sample correlation heatmap", "title")
    add_text(parts, 40, 60, "Pearson correlation based on filtered median-normalized log2 matrix; darker cells indicate higher correlation.", "subtitle")
    x0, y0 = 80, 105
    cell = 4.0
    size = cell * n
    values = corr.to_numpy()
    for i in range(n):
        for j in range(n):
            color = color_gradient(float(values[i, j]), 0.85, 1.0)
            parts.append(
                f'<rect x="{x0 + j * cell:.2f}" y="{y0 + i * cell:.2f}" '
                f'width="{cell:.2f}" height="{cell:.2f}" fill="{color}" shape-rendering="crispEdges"/>'
            )

    meta = sample_qc.set_index("sample_id")
    for j, sample in enumerate(ordered):
        color = HEATMAP_REGION_COLORS.get(meta.loc[sample, "region"], HEATMAP_REGION_COLORS["UNKNOWN"])
        parts.append(
            f'<rect x="{x0 + j * cell:.2f}" y="{y0 - 16:.2f}" width="{cell:.2f}" height="12" fill="{color}" shape-rendering="crispEdges"/>'
        )
        parts.append(
            f'<rect x="{x0 - 16:.2f}" y="{y0 + j * cell:.2f}" width="12" height="{cell:.2f}" fill="{color}" shape-rendering="crispEdges"/>'
        )
    parts.append(f'<rect x="{x0}" y="{y0}" width="{size}" height="{size}" fill="none" stroke="#333" stroke-width="1"/>')

    # Region boundaries.
    last_region = None
    for idx, sample in enumerate(ordered + [None]):
        region = meta.loc[sample, "region"] if sample is not None else None
        if idx == 0:
            last_region = region
        elif region != last_region:
            xpos = x0 + idx * cell
            ypos = y0 + idx * cell
            parts.append(f'<line x1="{xpos:.1f}" y1="{y0}" x2="{xpos:.1f}" y2="{y0+size}" stroke="#222" stroke-width="0.8"/>')
            parts.append(f'<line x1="{x0}" y1="{ypos:.1f}" x2="{x0+size}" y2="{ypos:.1f}" stroke="#222" stroke-width="0.8"/>')
            last_region = region

    # Color bar.
    cb_x, cb_y, cb_w, cb_h = 1370, 160, 28, 240
    for k in range(cb_h):
        val = 1.0 - k / cb_h * 0.15
        parts.append(f'<rect x="{cb_x}" y="{cb_y+k}" width="{cb_w}" height="1" fill="{color_gradient(val)}"/>')
    parts.append(f'<rect x="{cb_x}" y="{cb_y}" width="{cb_w}" height="{cb_h}" fill="none" stroke="#333"/>')
    add_text(parts, cb_x + 38, cb_y + 4, "1.00", "tick")
    add_text(parts, cb_x + 38, cb_y + cb_h, "0.85", "tick")
    add_text(parts, cb_x - 10, cb_y - 14, "r", "label")
    legend_x = 1280
    legend_regions = [r for r in HEATMAP_REGION_ORDER if r in set(sample_qc["region"])]
    for i, region in enumerate(legend_regions):
        yy = 460 + i * 18
        color = HEATMAP_REGION_COLORS.get(region, HEATMAP_REGION_COLORS["UNKNOWN"])
        parts.append(f'<rect x="{legend_x}" y="{yy-10}" width="11" height="11" fill="{color}"/>')
        add_text(parts, legend_x + 16, yy, region, "legend")
    save_svg(parts, out_path)


def figure_6_qc_cv(qc_cv: pd.DataFrame, out_path: Path) -> None:
    df = qc_cv.replace([np.inf, -np.inf], np.nan).dropna(subset=["qc_mean_log2_intensity", "qc_cv_percent"])
    df = df[df["qc_detected_fraction"] >= 0.8].copy()
    width, height = 1100, 820
    parts = svg_header(width, height)
    add_text(parts, 40, 38, "Figure 2E. QC technical variation by abundance", "title")
    add_text(parts, 40, 60, "Each point is a protein group detected in at least 80% of the 25 QC samples.", "subtitle")
    x0, y0, w, h = 105, 105, 780, 600
    xmin, xmax = float(df["qc_mean_log2_intensity"].min()), float(df["qc_mean_log2_intensity"].max())
    # Show the full QC acceptance range without clipping high-CV proteins at a percentile cap.
    ymax = 100.0
    ymin = 0.0
    add_axes(parts, x0, y0, w, h)
    add_y_ticks(parts, x0, y0, h, ymin, ymax, ticks=6, fmt="{:.0f}%")
    for i in range(6):
        value = xmin + (xmax - xmin) * i / 5
        x = scale(value, xmin, xmax, x0, x0 + w)
        parts.append(f'<line x1="{x:.1f}" y1="{y0+h}" x2="{x:.1f}" y2="{y0+h+5}" class="axis"/>')
        add_text(parts, x, y0 + h + 20, f"{value:.1f}", "tick", "middle")
    for threshold in [20, 30]:
        if threshold <= ymax:
            y = scale(threshold, ymin, ymax, y0 + h, y0)
            parts.append(
                f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+w}" y2="{y:.1f}" '
                'stroke="#b2182b" stroke-dasharray="5,5"/>'
            )
            add_text(parts, x0 + w + 8, y + 4, f"CV {threshold}%", "tick")

    sample = df
    for _, row in sample.iterrows():
        x = scale(row["qc_mean_log2_intensity"], xmin, xmax, x0, x0 + w)
        yval = min(row["qc_cv_percent"], ymax)
        y = scale(yval, ymin, ymax, y0 + h, y0)
        opacity = "0.30" if row["qc_cv_percent"] <= ymax else "0.12"
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.0" fill="#2c7fb8" opacity="{opacity}">'
            f'<title>{esc(row["PG.ProteinAccessions"])} CV={row["qc_cv_percent"]:.1f}%</title></circle>'
        )
    add_text(parts, x0 + w / 2, y0 + h + 50, "Mean log2 intensity in QC samples", "label", "middle")
    add_text(parts, 35, y0 + h / 2, "QC CV (%)", "label")

    total = len(df)
    below20 = int((df["qc_cv_percent"] <= 20).sum())
    below30 = int((df["qc_cv_percent"] <= 30).sum())
    add_text(parts, 915, 125, f"Proteins plotted: {total}", "legend")
    add_text(parts, 915, 150, f"CV <= 20%: {below20} ({below20 / total * 100:.1f}%)", "legend")
    add_text(parts, 915, 175, f"CV <= 30%: {below30} ({below30 / total * 100:.1f}%)", "legend")
    save_svg(parts, out_path)


def write_summary_md(
    out_path: Path,
    sample_qc: pd.DataFrame,
    qc_cv: pd.DataFrame,
    file_count: int,
) -> None:
    worst = sample_qc.sort_values("missing_percent", ascending=False).head(8)
    qc_only = sample_qc[sample_qc["sample_type"] == "QC"]
    cv_filtered = qc_cv[(qc_cv["qc_detected_fraction"] >= 0.8) & qc_cv["qc_cv_percent"].notna()]
    below20 = (cv_filtered["qc_cv_percent"] <= 20).mean() * 100
    below30 = (cv_filtered["qc_cv_percent"] <= 30).mean() * 100
    lines = [
        "# QC Analysis Summary",
        "",
        "This report was generated by `04_reproducible_scripts/qc_svg_analysis.py`.",
        "",
        "## Input Overview",
        "",
        f"- Files analyzed: {file_count}",
        f"- Samples analyzed: {sample_qc.shape[0]}",
        f"- QC samples: {qc_only.shape[0]}",
        "",
        "## QC Sample Stability",
        "",
        f"- QC median detected protein groups: {qc_only['detected_protein_groups'].median():.0f}",
        f"- QC median missing rate: {qc_only['missing_percent'].median():.2f}%",
        f"- QC proteins with CV <= 20%: {below20:.1f}% among proteins detected in at least 80% QC samples",
        f"- QC proteins with CV <= 30%: {below30:.1f}% among proteins detected in at least 80% QC samples",
        "",
        "## Highest Missing-Rate Samples",
        "",
        "| sample_id | region | group | missing_percent | detected_protein_groups |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in worst.iterrows():
        lines.append(
            f"| {row['sample_id']} | {row['region']} | {row['group']} | "
            f"{row['missing_percent']:.2f} | {int(row['detected_protein_groups'])} |"
        )
    lines.extend(
        [
            "",
            "## Generated Figures",
            "",
            "- `figures/figure_2D_sample_detection_missing.svg`",
            "- `figures/figure_2E_qc_cv_vs_abundance.svg`",
            "- `figures/figure_2F_sample_correlation_heatmap.svg`",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate the three selected QC figures for DIA-LFQ reports.")
    parser.add_argument("--data-dir", default="data", help="Directory containing CSV report files.")
    parser.add_argument("--out-dir", default="results/qc", help="Output directory.")
    parser.add_argument("--min-detection-fraction", type=float, default=0.8, help="Minimum detection fraction for the correlation matrix.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    table_dir = out_dir / "tables"
    safe_mkdir(fig_dir)
    safe_mkdir(table_dir)
    for obsolete_name in [
        "figure_2_qc_injection_stability.svg",
        "figure_3_intensity_distribution_normalization.svg",
        "figure_4_pca.svg",
        "figure_1_sample_detection_missing.svg",
        "figure_5_sample_correlation_heatmap.svg",
        "figure_6_qc_cv_vs_abundance.svg",
    ]:
        (fig_dir / obsolete_name).unlink(missing_ok=True)
    for obsolete_name in [
        "protein_detection_summary_by_file.csv",
        "qc_injection_stability_metrics.csv",
        "pca_scores.csv",
    ]:
        (table_dir / obsolete_name).unlink(missing_ok=True)

    reports = load_all_reports(data_dir)
    matrix, sample_meta, annotation = build_combined_matrix(reports)
    sample_qc = sample_qc_summary_from_reports(reports)
    sample_meta = add_analysis_decisions(sample_meta)
    sample_qc = add_analysis_decisions(sample_qc)
    qc_cv = qc_cv_table(matrix, sample_meta, annotation)
    complete = prepare_complete_matrix(matrix, args.min_detection_fraction)
    corr = pearson_corr_matrix(complete)

    sample_meta.to_csv(table_dir / "sample_metadata.csv", index=False)
    sample_qc.to_csv(table_dir / "sample_qc_summary.csv", index=False)
    qc_cv.to_csv(table_dir / "qc_cv_by_protein.csv", index=False)
    corr.to_csv(table_dir / "sample_correlation_matrix.csv")

    figure_1_sample_overview(sample_qc, fig_dir / "figure_2D_sample_detection_missing.svg")
    figure_6_qc_cv(qc_cv, fig_dir / "figure_2E_qc_cv_vs_abundance.svg")
    figure_5_correlation_heatmap(corr, sample_qc, fig_dir / "figure_2F_sample_correlation_heatmap.svg")
    write_summary_md(out_dir / "qc_analysis_summary.md", sample_qc, qc_cv, len(reports))

    print(f"Loaded reports: {len(reports)}")
    print(f"Combined matrix: {matrix.shape[0]} proteins x {matrix.shape[1]} samples")
    print(f"Filtered matrix for correlation: {complete.shape[0]} proteins x {complete.shape[1]} samples")
    print(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
