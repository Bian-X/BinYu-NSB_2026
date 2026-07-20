#!/usr/bin/env python
"""Create the three selected standalone QC PDF panels.

The report uses only reportlab plus pandas/numpy so it works in the bundled
runtime without matplotlib, seaborn, sklearn, or SVG conversion libraries.
The main entry point is `run_qc_pdf.py`.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


REGION_ORDER = ["QC", "PL", "IL", "SSC", "MSC", "BLA", "CEA", "DHPC", "PVH", "VHPC", "PMD", "LPB"]
REGION_COLORS = {
    "QC": colors.HexColor("#2f3437"),
    "BLA": colors.HexColor("#6a3d9a"),
    "CEA": colors.HexColor("#c2a5cf"),
    "IL": colors.HexColor("#1f78b4"),
    "PL": colors.HexColor("#a6cee3"),
    "DHPC": colors.HexColor("#1b7837"),
    "VHPC": colors.HexColor("#a6dba0"),
    "SSC": colors.HexColor("#e66101"),
    "MSC": colors.HexColor("#fdb863"),
    "PVH": colors.HexColor("#d73027"),
    "PMD": colors.HexColor("#8c510a"),
    "LPB": colors.HexColor("#636363"),
    "UNKNOWN": colors.HexColor("#8c8c8c"),
}

HEATMAP_REGION_ORDER = ["QC", "BLA", "CEA", "IL", "PL", "DHPC", "VHPC", "SSC", "MSC", "PVH", "PMD", "LPB"]
# Heatmap annotations use the global brain-region palette used by all QC figures.
HEATMAP_REGION_COLORS = REGION_COLORS


def scale(value: float, old_min: float, old_max: float, new_min: float, new_max: float) -> float:
    if not np.isfinite(value) or old_max == old_min:
        return (new_min + new_max) / 2
    return new_min + (value - old_min) * (new_max - new_min) / (old_max - old_min)


def region_color(region: str):
    return REGION_COLORS.get(str(region), REGION_COLORS["UNKNOWN"])


def draw_title(c: canvas.Canvas, title: str, subtitle: str, page_num: int, width: float, height: float) -> None:
    c.setFillColor(colors.HexColor("#222222"))
    c.setFont("Helvetica-Bold", 15)
    c.drawString(0.45 * inch, height - 0.45 * inch, title)
    c.setFont("Helvetica", 8.5)
    c.setFillColor(colors.HexColor("#555555"))
    c.drawString(0.45 * inch, height - 0.63 * inch, subtitle)


def draw_axes(c: canvas.Canvas, x: float, y: float, w: float, h: float) -> None:
    c.setStrokeColor(colors.HexColor("#333333"))
    c.setLineWidth(0.6)
    c.line(x, y, x + w, y)
    c.line(x, y, x, y + h)


def draw_y_ticks(
    c: canvas.Canvas,
    x: float,
    y: float,
    w: float,
    h: float,
    vmin: float,
    vmax: float,
    n: int = 5,
    fmt: str = "{:.0f}",
) -> None:
    c.setFont("Helvetica", 6.5)
    c.setStrokeColor(colors.HexColor("#dddddd"))
    for i in range(n):
        value = vmin + (vmax - vmin) * i / (n - 1)
        yy = scale(value, vmin, vmax, y, y + h)
        c.setStrokeColor(colors.HexColor("#dddddd"))
        c.line(x, yy, x + w, yy)
        c.setStrokeColor(colors.HexColor("#333333"))
        c.line(x - 3, yy, x, yy)
        c.setFillColor(colors.HexColor("#555555"))
        c.drawRightString(x - 5, yy - 2, fmt.format(value))


def draw_legend(c: canvas.Canvas, x: float, y: float, regions: list[str]) -> None:
    c.setFont("Helvetica", 7)
    for i, region in enumerate(regions):
        yy = y - i * 12
        c.setFillColor(region_color(region))
        c.rect(x, yy - 7, 7, 7, fill=1, stroke=0)
        c.setFillColor(colors.HexColor("#333333"))
        c.drawString(x + 10, yy - 6, region)


def page_1_sample_overview(c: canvas.Canvas, sample_qc: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 2D. Sample detection depth and missingness overview",
        "Bars show detected protein groups; points show sample missing percentage.",
        1,
        width,
        height,
    )
    df = sample_qc.sort_values(["region_rank", "group_rank", "replicate", "filename_index", "qc_index"]).reset_index(drop=True)
    x0, y0, w, h = 0.7 * inch, 4.55 * inch, 9.0 * inch, 2.0 * inch
    y1, h2 = 1.25 * inch, 2.0 * inch
    draw_axes(c, x0, y0, w, h)
    det_min, det_max = 0, float(df["detected_protein_groups"].max() * 1.05)
    draw_y_ticks(c, x0, y0, w, h, det_min, det_max, 5)
    n = len(df)
    bar_w = max(0.35, w / n * 0.82)
    region_boundaries = [
        index for index in range(1, n) if df.loc[index, "region"] != df.loc[index - 1, "region"]
    ]
    c.setDash(3, 3)
    c.setStrokeColor(colors.HexColor("#9e9e9e"))
    c.setLineWidth(0.5)
    for index in region_boundaries:
        x = x0 + index * w / n
        c.line(x, y1, x, y0 + h)
    c.setDash()
    for i, row in df.iterrows():
        cx = x0 + i * w / n
        bar_h = scale(row["detected_protein_groups"], det_min, det_max, 0, h)
        c.setFillColor(region_color(row["region"]))
        c.rect(cx, y0, bar_w, bar_h, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#333333"))
    c.setFont("Helvetica", 7)
    c.drawString(0.28 * inch, y0 + h / 2, "Detected")

    draw_axes(c, x0, y1, w, h2)
    miss_min, miss_max = 0, max(10.0, float(df["missing_percent"].max() * 1.1))
    draw_y_ticks(c, x0, y1, w, h2, miss_min, miss_max, 6, "{:.0f}%")
    for i, row in df.iterrows():
        cx = x0 + i * w / n + bar_w / 2
        cy = scale(row["missing_percent"], miss_min, miss_max, y1, y1 + h2)
        c.setFillColor(region_color(row["region"]))
        radius = 2.3 if row["sample_type"] == "QC" else 1.6
        c.circle(cx, cy, radius, fill=1, stroke=0)
    threshold_y = scale(10, miss_min, miss_max, y1, y1 + h2)
    c.setDash(4, 3)
    c.setStrokeColor(colors.HexColor("#b2182b"))
    c.line(x0, threshold_y, x0 + w, threshold_y)
    c.setDash()
    c.setFillColor(colors.HexColor("#b2182b"))
    c.setFont("Helvetica", 6.5)
    c.drawString(x0 + w + 4, threshold_y - 2, "10%")
    c.setFillColor(colors.HexColor("#333333"))
    c.drawString(0.28 * inch, y1 + h2 / 2, "Missing")
    draw_legend(c, 9.95 * inch, 6.25 * inch, [r for r in REGION_ORDER if r in set(df["region"])])


def page_2_qc_injection(c: canvas.Canvas, qc_metrics: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 2. QC injection-order stability",
        "The 25 QC samples are ordered by QC injection index.",
        2,
        width,
        height,
    )
    df = qc_metrics.sort_values("qc_index").reset_index(drop=True)
    panels = [
        ("detected_protein_groups", "Detected proteins", "{:.0f}"),
        ("missing_percent", "Missing rate (%)", "{:.1f}"),
        ("median_log2_intensity", "Median log2 intensity", "{:.1f}"),
        ("correlation_to_qc_median", "Correlation to QC median", "{:.3f}"),
    ]
    positions = [(0.8 * inch, 4.35 * inch), (5.9 * inch, 4.35 * inch), (0.8 * inch, 1.1 * inch), (5.9 * inch, 1.1 * inch)]
    panel_w, panel_h = 4.35 * inch, 2.25 * inch
    for (col, title, fmt), (x0, y0) in zip(panels, positions):
        c.setFillColor(colors.HexColor("#333333"))
        c.setFont("Helvetica-Bold", 8)
        c.drawString(x0, y0 + panel_h + 11, title)
        draw_axes(c, x0, y0, panel_w, panel_h)
        values = df[col].astype(float)
        vmin, vmax = float(values.min()), float(values.max())
        pad = (vmax - vmin) * 0.15 if vmax > vmin else 1
        vmin = 0 if col == "missing_percent" else vmin - pad
        vmax += pad
        draw_y_ticks(c, x0, y0, panel_w, panel_h, vmin, vmax, 5, fmt)
        points = []
        for _, row in df.iterrows():
            x = scale(row["qc_index"], 1, 25, x0, x0 + panel_w)
            y = scale(row[col], vmin, vmax, y0, y0 + panel_h)
            points.append((x, y))
        c.setStrokeColor(colors.HexColor("#2f3437"))
        c.setLineWidth(1.2)
        for (x1, y1), (x2, y2) in zip(points, points[1:]):
            c.line(x1, y1, x2, y2)
        c.setFillColor(colors.HexColor("#2f3437"))
        for x, y in points:
            c.circle(x, y, 2.2, fill=1, stroke=0)
        c.setFont("Helvetica", 6.5)
        for tick in [1, 5, 10, 15, 20, 25]:
            x = scale(tick, 1, 25, x0, x0 + panel_w)
            c.drawCentredString(x, y0 - 12, str(tick))


def page_3_intensity(c: canvas.Canvas, sample_qc: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 3. Log2 intensity distribution before and after normalization",
        "Each vertical segment is one sample IQR; the dot is the sample median.",
        3,
        width,
        height,
    )
    df = sample_qc.sort_values(["region_rank", "group_rank", "replicate", "filename_index", "qc_index"]).reset_index(drop=True)
    # After median normalization, all sample medians shift to the global median.
    global_median = df["median_log2_intensity"].median()
    after = df.copy()
    shifts = global_median - after["median_log2_intensity"]
    after["q25_log2_intensity"] = after["q25_log2_intensity"] + shifts
    after["q75_log2_intensity"] = after["q75_log2_intensity"] + shifts
    after["median_log2_intensity"] = global_median
    combined_min = min(df["q25_log2_intensity"].min(), after["q25_log2_intensity"].min()) - 0.5
    combined_max = max(df["q75_log2_intensity"].max(), after["q75_log2_intensity"].max()) + 0.5
    for title, data, y0 in [
        ("Before median normalization", df, 4.25 * inch),
        ("After median normalization", after, 1.25 * inch),
    ]:
        x0, w, h = 0.75 * inch, 9.0 * inch, 2.0 * inch
        c.setFillColor(colors.HexColor("#333333"))
        c.setFont("Helvetica-Bold", 8)
        c.drawString(x0, y0 + h + 11, title)
        draw_axes(c, x0, y0, w, h)
        draw_y_ticks(c, x0, y0, w, h, combined_min, combined_max, 5, "{:.1f}")
        n = len(data)
        for i, row in data.iterrows():
            x = x0 + i * w / n
            y25 = scale(row["q25_log2_intensity"], combined_min, combined_max, y0, y0 + h)
            y75 = scale(row["q75_log2_intensity"], combined_min, combined_max, y0, y0 + h)
            ym = scale(row["median_log2_intensity"], combined_min, combined_max, y0, y0 + h)
            c.setStrokeColor(region_color(row["region"]))
            c.setLineWidth(0.45)
            c.line(x, y25, x, y75)
            c.setFillColor(region_color(row["region"]))
            c.circle(x, ym, 1.4 if row["sample_type"] != "QC" else 2.0, fill=1, stroke=0)
    draw_legend(c, 9.95 * inch, 6.2 * inch, [r for r in REGION_ORDER if r in set(df["region"])])


def page_4_pca(c: canvas.Canvas, pca: pd.DataFrame, sample_qc: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 4. PCA of biological samples and QC samples",
        "PCA uses median-normalized log2 protein-group quantities with row-median imputation.",
        4,
        width,
        height,
    )
    df = sample_qc.merge(pca, on="sample_id", how="inner")
    x0, y0, w, h = 0.95 * inch, 1.0 * inch, 7.4 * inch, 5.7 * inch
    xmin, xmax = df["PC1"].min(), df["PC1"].max()
    ymin, ymax = df["PC2"].min(), df["PC2"].max()
    xpad, ypad = (xmax - xmin) * 0.08, (ymax - ymin) * 0.08
    xmin, xmax, ymin, ymax = xmin - xpad, xmax + xpad, ymin - ypad, ymax + ypad
    draw_axes(c, x0, y0, w, h)
    draw_y_ticks(c, x0, y0, w, h, ymin, ymax, 6, "{:.1f}")
    c.setFont("Helvetica", 6.5)
    for i in range(6):
        val = xmin + (xmax - xmin) * i / 5
        x = scale(val, xmin, xmax, x0, x0 + w)
        c.drawCentredString(x, y0 - 12, f"{val:.1f}")
    for _, row in df.iterrows():
        x = scale(row["PC1"], xmin, xmax, x0, x0 + w)
        y = scale(row["PC2"], ymin, ymax, y0, y0 + h)
        c.setFillColor(region_color(row["region"]))
        if row["sample_type"] == "QC":
            c.rect(x - 2.5, y - 2.5, 5, 5, fill=1, stroke=0)
        else:
            c.circle(x, y, 2.1, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#333333"))
    c.setFont("Helvetica", 8)
    c.drawCentredString(x0 + w / 2, y0 - 31, "PC1")
    c.drawString(0.35 * inch, y0 + h / 2, "PC2")
    draw_legend(c, 8.65 * inch, 6.3 * inch, [r for r in REGION_ORDER if r in set(df["region"])])


def corr_color(value: float):
    vmin, vmax = 0.85, 1.0
    value = max(vmin, min(vmax, float(value)))
    t = (value - vmin) / (vmax - vmin)
    r = (222 * (1 - t) + 8 * t) / 255
    g = (235 * (1 - t) + 69 * t) / 255
    b = (247 * (1 - t) + 148 * t) / 255
    return colors.Color(r, g, b)


def page_5_corr(c: canvas.Canvas, corr: pd.DataFrame, sample_qc: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 2F. Sample correlation heatmap",
        "Pearson correlation based on filtered median-normalized log2 matrix; darker cells indicate higher correlation.",
        5,
        width,
        height,
    )
    heatmap_region_rank = {region: index for index, region in enumerate(HEATMAP_REGION_ORDER)}
    ordered = sample_qc.assign(
        heatmap_region_rank=sample_qc["region"].map(heatmap_region_rank).fillna(len(HEATMAP_REGION_ORDER))
    ).sort_values(["heatmap_region_rank", "group_rank", "replicate", "filename_index", "qc_index"])["sample_id"].tolist()
    ordered = [s for s in ordered if s in corr.index]
    corr = corr.loc[ordered, ordered]
    n = len(ordered)
    x0, y0 = 0.75 * inch, 0.45 * inch
    cell = 1.8
    size = cell * n
    vals = corr.to_numpy()
    c.setLineWidth(0)
    for i in range(n):
        yy = y0 + size - (i + 1) * cell
        for j in range(n):
            c.setFillColor(corr_color(vals[i, j]))
            c.rect(x0 + j * cell, yy, cell, cell, fill=1, stroke=0)
    meta = sample_qc.set_index("sample_id")
    for j, sample in enumerate(ordered):
        c.setFillColor(HEATMAP_REGION_COLORS.get(meta.loc[sample, "region"], HEATMAP_REGION_COLORS["UNKNOWN"]))
        c.rect(x0 + j * cell, y0 + size + 3, cell, 6, fill=1, stroke=0)
        c.rect(x0 - 9, y0 + size - (j + 1) * cell, 6, cell, fill=1, stroke=0)
    c.setStrokeColor(colors.HexColor("#333333"))
    c.setLineWidth(0.6)
    c.rect(x0, y0, size, size, fill=0, stroke=1)
    # Color bar.
    cb_x, cb_y, cb_w, cb_h = 8.1 * inch, 4.75 * inch, 0.22 * inch, 1.65 * inch
    steps = 80
    for k in range(steps):
        val = 1.0 - k / steps * 0.15
        c.setFillColor(corr_color(val))
        c.rect(cb_x, cb_y + cb_h * (steps - k - 1) / steps, cb_w, cb_h / steps + 0.2, fill=1, stroke=0)
    c.setStrokeColor(colors.HexColor("#333333"))
    c.rect(cb_x, cb_y, cb_w, cb_h, fill=0, stroke=1)
    c.setFont("Helvetica", 7)
    c.setFillColor(colors.HexColor("#333333"))
    c.drawString(cb_x + cb_w + 5, cb_y + cb_h - 3, "1.00")
    c.drawString(cb_x + cb_w + 5, cb_y - 2, "0.85")
    c.drawString(cb_x, cb_y + cb_h + 10, "r")
    c.setFont("Helvetica", 7)
    for i, region in enumerate([r for r in HEATMAP_REGION_ORDER if r in set(sample_qc["region"])]):
        yy = 4.25 * inch - i * 12
        c.setFillColor(HEATMAP_REGION_COLORS.get(region, HEATMAP_REGION_COLORS["UNKNOWN"]))
        c.rect(8.1 * inch, yy - 7, 7, 7, fill=1, stroke=0)
        c.setFillColor(colors.HexColor("#333333"))
        c.drawString(8.1 * inch + 10, yy - 6, region)


def page_6_qc_cv(c: canvas.Canvas, qc_cv: pd.DataFrame, width: float, height: float) -> None:
    draw_title(
        c,
        "Figure 2E. QC technical variation by abundance",
        "Each point is a protein group detected in at least 80% of the 25 QC samples.",
        6,
        width,
        height,
    )
    df = qc_cv.replace([np.inf, -np.inf], np.nan).dropna(subset=["qc_mean_log2_intensity", "qc_cv_percent"])
    df = df[df["qc_detected_fraction"] >= 0.8].copy()
    x0, y0, w, h = 0.9 * inch, 1.0 * inch, 7.2 * inch, 5.55 * inch
    xmin, xmax = df["qc_mean_log2_intensity"].min(), df["qc_mean_log2_intensity"].max()
    ymax = 100.0
    ymin = 0
    draw_axes(c, x0, y0, w, h)
    draw_y_ticks(c, x0, y0, w, h, ymin, ymax, 6, "{:.0f}%")
    c.setFont("Helvetica", 6.5)
    for i in range(6):
        val = xmin + (xmax - xmin) * i / 5
        x = scale(val, xmin, xmax, x0, x0 + w)
        c.drawCentredString(x, y0 - 12, f"{val:.1f}")
    for threshold in [20, 30]:
        if threshold <= ymax:
            yy = scale(threshold, ymin, ymax, y0, y0 + h)
            c.setDash(4, 3)
            c.setStrokeColor(colors.HexColor("#b2182b"))
            c.line(x0, yy, x0 + w, yy)
            c.setDash()
            c.setFillColor(colors.HexColor("#b2182b"))
            c.drawString(x0 + w + 4, yy - 2, f"CV {threshold}%")
    c.setFillColor(colors.Color(0.17, 0.5, 0.72, alpha=0.22))
    for _, row in df.iterrows():
        x = scale(row["qc_mean_log2_intensity"], xmin, xmax, x0, x0 + w)
        y = scale(min(row["qc_cv_percent"], ymax), ymin, ymax, y0, y0 + h)
        c.circle(x, y, 1.2, fill=1, stroke=0)
    total = len(df)
    below20 = int((df["qc_cv_percent"] <= 20).sum())
    below30 = int((df["qc_cv_percent"] <= 30).sum())
    c.setFillColor(colors.HexColor("#333333"))
    c.setFont("Helvetica", 8)
    c.drawString(8.35 * inch, 6.05 * inch, f"Proteins plotted: {total}")
    c.drawString(8.35 * inch, 5.82 * inch, f"CV <= 20%: {below20} ({below20 / total * 100:.1f}%)")
    c.drawString(8.35 * inch, 5.59 * inch, f"CV <= 30%: {below30} ({below30 / total * 100:.1f}%)")
    c.drawCentredString(x0 + w / 2, y0 - 31, "Mean log2 intensity in QC samples")
    c.drawString(0.35 * inch, y0 + h / 2, "QC CV (%)")


def _write_single_page_pdf(out_pdf: Path, draw_page) -> None:
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    page_size = landscape(letter)
    width, height = page_size
    c = canvas.Canvas(str(out_pdf), pagesize=page_size)
    draw_page(c, width, height)
    c.save()


def create_figure_pdfs(
    sample_qc: pd.DataFrame,
    qc_cv: pd.DataFrame,
    corr: pd.DataFrame,
    out_dir: Path,
) -> list[Path]:
    """Create the three selected Figure 2 QC panels as separate PDFs."""
    outputs = [
        out_dir / "figure_2D_detection_missing.pdf",
        out_dir / "figure_2E_qc_cv_abundance.pdf",
        out_dir / "figure_2F_sample_correlation.pdf",
    ]
    _write_single_page_pdf(outputs[0], lambda c, w, h: page_1_sample_overview(c, sample_qc, w, h))
    _write_single_page_pdf(outputs[1], lambda c, w, h: page_6_qc_cv(c, qc_cv, w, h))
    _write_single_page_pdf(outputs[2], lambda c, w, h: page_5_corr(c, corr, sample_qc, w, h))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Create three standalone Figure 2 QC PDFs from QC tables.")
    parser.add_argument("--table-dir", default="results/qc/tables", help="QC table directory.")
    parser.add_argument("--out-dir", default="output/pdf", help="Directory for the three PDF panels.")
    args = parser.parse_args()
    table_dir = Path(args.table_dir)
    outputs = create_figure_pdfs(
        pd.read_csv(table_dir / "sample_qc_summary.csv"),
        pd.read_csv(table_dir / "qc_cv_by_protein.csv"),
        pd.read_csv(table_dir / "sample_correlation_matrix.csv", index_col=0),
        Path(args.out_dir),
    )
    for output in outputs:
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
