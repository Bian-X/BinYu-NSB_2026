# Input data directory

This directory is intentionally empty in the public release.

To reproduce the Figure 2 QC panels, add the twelve original Spectronaut DIA-LFQ protein-group CSV reports here: one `..._QC_DIA_...csv` report plus eleven regional `..._<REGION>_DIA_...csv` reports. The reports are read directly; no file needs to be renamed.

Required columns in every report:

- `PG.ProteinAccessions`
- one or more columns ending in `.PG.Quantity`

Do not commit the reports, generated supporting tables, or figures unless their public release has been approved.

