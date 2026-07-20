# Figure 2 QC manuscript-alignment record

The supported PDF workflow was rerun against the archived 12 input reports on 2026-07-20. The results below match the revised manuscript's Figure 2 legend and Sections 3.10.2-3.10.4.

| Manuscript statement | Recomputed result | Status |
|---|---:|---|
| Measurements / protein groups in combined matrix | 289 / 11,322 | Match |
| Figure 2D detected protein groups (range; median [IQR]) | 6,143-9,294; 8,940 [8,852-9,005] | Match |
| Figure 2D missingness (median [IQR]) | 2.61% [2.05%-3.44%] | Match |
| Extreme low-coverage PVH sample | 6,143 protein groups; 33.21% missing | Match |
| Retained CEA sample | 10.87% missing | Match |
| QC report rows / injections | 9,093 / 25 | Match |
| QC detected protein groups (range; median [IQR]) | 8,924-9,013; 8,980 [8,963-8,985] | Match |
| QC missingness (range; median [IQR]) | 0.88%-1.86%; 1.24% [1.19%-1.43%] | Match |
| Figure 2E proteins retained at >=80% QC detection | 8,908 | Match |
| CV <=20% / <=30% | 6,957 (78.1%) / 8,050 (90.4%) | Match |
| Proteins retained for Figure 2F correlation | 8,072 | Match |
| Median BLA-CEA / IL-PL / dHPC-vHPC / sSC-mSC correlations | 0.9572 / 0.9583 / 0.9617 / 0.9550 | Match |

## Important scope note

The QC code records the PVH sample `Astral_MSP2600441_253_G6_1_PVH` as excluded from downstream biological analyses, while Figure 2D and Figure 2F deliberately retain all 289 acquisitions as described in the manuscript. This is consistent: those panels are acquisition-level QC descriptions, not the downstream differential-analysis input set.

## Documentation correction made for the public package

The archived QC folder included both a current three-PDF entry point and an older SVG-generating command with outdated paths in one summary file. The public package designates the PDF entry point as the sole supported workflow and excludes obsolete output-path references from its documentation. The archival SVG routine remains in `scripts/qc_analysis.py` but is not used by the release instructions.

