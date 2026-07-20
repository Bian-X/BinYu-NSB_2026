# GitHub release checklist

Before making this repository public:

1. Confirm that `data/` contains no input CSVs, raw MS files, restricted sample
   metadata, credentials, or participant/animal identifiers beyond what is
   approved for release.
2. Run the workflow with a permitted copy of the input data and retain the
   resulting `environment/sessionInfo.txt` as the execution record.
3. Review `git status --ignored` and remove generated `figures/` and `results/`
   files unless you intentionally want to distribute them as versioned expected
   outputs.
4. Choose and add a license. Do not assume an open-source license without the
   agreement of the copyright holders and data owners.
5. Create a GitHub release or archive a tagged version with Zenodo, then add the
   URL or DOI to the manuscript's Code availability statement.
6. Replace any manuscript placeholders for repository URL, archival DOI, access
   restrictions, and software license.

Local validation without input data:

```bash
python -m py_compile scripts/python/brain_region_pca.py
bash -n run_all.sh
```

With R available, additionally parse the R scripts:

```bash
Rscript -e "invisible(lapply(list.files('scripts/R', full.names=TRUE), parse))"
```
