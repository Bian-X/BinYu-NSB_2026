# Analysis methods

## Inputs and identity

Eleven regional processed protein-group reports contain 24 biological measurements
each (G1–G6, four mice per group). The pooled-QC report contributes 25 acquisitions.
Animal identity is shared across regions; 264 regional measurements are not 264
independent mice. G1: untreated; G2: unpaired tone/shock; G3: US-only; G4: paired;
G5: contextual retrieval; G6: tone-cued retrieval. Inspect the supplied sample
metadata and the original study design before interpreting these conditions.

## Primary preprocessing and differential abundance

Each region is processed separately. Retain protein groups whose `PG.Organisms`
contains `Mus musculus`, exclude the exact semicolon-delimited accession
`iRT-Kit_WR_fusion`, and treat nonpositive or nonfinite quantities as missing.
Exclude the low-coverage PVH/G6 measurement in `sample_analysis/sample_exclusions.csv` before
filtering. Require at least 20 observed measurements (20/24, or 20/23 for primary
PVH). Log2-transform, normalize by observed sample medians, and fill remaining
missing cells with the protein-group median on that normalized log2 scale.
Normalization subtracts each observed sample median and adds the median of
those sample medians. Completion is deterministic and occurs after filtering
and normalization.

Fit all six groups with limma (`trend=TRUE`, `robust=TRUE`). Extract G5−G1 and,
for the selected GO example, BLA G4−G1. Candidates require BH P<0.05 within the
region/contrast and |log2FC|>log2(1.5). Output observed counts support interpretation.
This does not establish FDR control across all regions and all contrasts.

The primary analysis has 263 regional measurements from 24 mice. The include-all
264-measurement sensitivity workflow repeats preprocessing from the input reports.
The excluded PVH/G6 measurement has 33.21% within-report missingness.
The 15% missingness exclusion rule was retrospective, not preregistered.

## GO over-representation

Map protein groups to a single unique Entrez gene for the primary analysis;
apply the same mapping to foreground and background, then deduplicate genes.
Foreground comes from the corresponding differential candidates, background from
tested protein groups. GOALL annotations define eligible BP/CC/MF terms with
5–500 background genes. Use upper-tail hypergeometric P values and one pooled
BH family containing all eligible terms, including zero-overlap terms (P=1).
Expanded multi-gene mapping is explicitly a sensitivity analysis.

PL G5−G1: 8,669 tested groups; 251 candidates (153 higher in G5, 98 lower);
242 foreground/8,449 background unique-mapped genes; 9,135 GO tests; zero BH hits.
Expanded PL mapping has 17 BH hits and is sensitive to ambiguous protein groups.
BLA G4−G1 unique mapping: 400 foreground/8,368 background genes; 9,097 terms;
37 within-family BH hits. Its foreground is not the PL volcano candidate set.

The GO example was selected during an exploratory screen of 11 regions × 15
contrasts. A supplementary BH calculation across all 1,510,215 screen GO tests
retains 23 BLA G4−G1 terms. Complete-case refitting gives 41 terms and retains
the original 37. These are sensitivity diagnostics, not prospective confirmation.
Run the optional screen and complete-case refit with `python run_all.py
--data-dir /path/to/processeddata --selection-sensitivity --rscript Rscript`.
The scripts are `compare_workflow_examples.R` and `validate_workflow_shortlist.R`
under sample_analysis/scripts/R/.

## PCA, heatmap and regional correlations

PCA inner-joins exact protein-group accessions across regions (7,390 groups),
centers features, does not scale to unit variance, and uses full SVD. Main PC1/PC2
explain 37.3211%/15.5856%; all measurements are shown as small circles.
The plotting area and Figure 3C page are square, with equal PC1/PC2 distance
units; the brain-region legend is outside the plotting area on the right.
Each regional ellipse uses the sample covariance and the 95% chi-square
contour with two degrees of freedom (scale sqrt(5.991464547)). This is an
approximate 95% data ellipse, not a confidence interval for the regional mean.

The G4−G2 heatmap uses an independent two-group model, at least 3/4 observations
per group, log2 protein-row median completion, `trend=TRUE`, `robust=FALSE`, and
nominal-P selection by default. BH values are always exported; use
`--selection-pvalue adjusted` with the standalone heatmap script for BH selection.
The exported `Status` column follows the selected P-value criterion.
This two-group branch does not apply sample-median normalization.
Missing region–group means are completed with the protein-row median before
row Z-scores. VisualScore-based
display selection is exploratory and does not imply FDR-significant pathways.

G5−G1 regional Pearson correlations use the union of regional candidates,
pairwise-complete values with at least 30 shared groups, and complete linkage on
1−r. Selection and coverage affect these correlations; they are not neural
connectivity or paired-animal hypothesis tests. Neither contrast isolates a pure
retrieval effect by itself.

## Acquisition QC and missingness

Retain all 289 acquisitions, including the excluded biological PVH measurement.
Detection, missingness, signal and CV use observed values without normalization
or completion. Within-report missingness uses the protein-group rows in that
report as its denominator. CV is 100 times the sample standard deviation divided
by the original-scale mean and requires at least 20 observed QC acquisitions.
Correlation matrices use one shared matrix:
at least 80% observed, log2, observed-sample median normalization, and protein-row
median completion. Region/injection/QC views reorder this same completed matrix.
Temporal slopes and R² are descriptive because common-reference and pairwise
correlations are dependent observations. Absence of a trend does not establish
absence of drift or confounding. Original event and retention-time records are
not synthesized when unavailable. Figure 2F displays the same matrix ordered by
region, group, and replicate; the standalone injection-order view uses metadata
order. Four CV values above 100% are drawn at the Figure 2E axis limit without
altering their exported values.

## QC populations

CV eligibility applies to all exported entries with >=20/25 observed pooled-QC
values (8,908). Injection-order regression additionally requires PG.Organisms
to contain Mus musculus (8,880). The 28 excluded entries consist of one iRT
standard and 27 entries annotated Unknown; organism assignment is not inferred
from free-text descriptions. The row-level audit is qc_regression_selection_audit.csv.
Regression uses observed log2 values, each row's available observations, and
the full QC injection-order span of 288 units. Fitted high/low ratios are
2**(abs(slope)*span); these processed quantities do not establish absence of
nonlinear drift, interruptions or condition-order confounding.
