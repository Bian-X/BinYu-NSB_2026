# Input data directory

This directory is intentionally empty in the public code repository.

Copy the 11 non-QC region-level CSV files here before running the workflow, or
provide their location with `-DataDir` (PowerShell) or `DATA_DIR` (Bash).

Each file must have a name matching the Spectronaut-style pattern
`_*REGION*_DIA_LFQ*.csv`, where `REGION` is one of:

`PMD`, `LPB`, `VHPC`, `PVH`, `BLA`, `DHPC`, `CEA`, `MSC`, `SSC`, `IL`, `PL`.

Required annotation columns:

- `PG.ProteinAccessions`
- `PG.Genes`
- `PG.Organisms`
- `PG.ProteinDescriptions`
- `PG.NrOfStrippedSequencesIdentified (Experiment-wide)`

The differential analyses require four columns matching each of the G2 and G4
replicates, with names ending in `.raw.PG.Quantity`. The PCA workflow expects
24 quantitative columns in each region report.

Do not commit input data, raw instrument files, or sample metadata unless their
distribution has been approved by the data owner and ethics/data-access terms.
