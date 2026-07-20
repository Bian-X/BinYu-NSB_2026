#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_executable="${PYTHON:-python3}"
rscript_executable="${RSCRIPT:-Rscript}"
font_family="${PROTEOMICS_FIGURE_FONT:-Arial}"
data_dir="${DATA_DIR:-}"
figure_dir="${root}/figures"
result_dir="${root}/results"

has_region_data() {
  local path="$1"
  [[ -d "${path}" ]] || return 1
  find "${path}" -maxdepth 1 -type f -name "*.csv" | grep -Eq "_(PMD|LPB|VHPC|PVH|BLA|DHPC|CEA|MSC|SSC|IL|PL)_DIA_LFQ.*\\.csv$"
}

if [[ -z "${data_dir}" ]]; then
  for candidate in "${root}/data" "${root}/data/data1" "${root}/../data/data1" "${root}/../data"; do
    if has_region_data "${candidate}"; then
      data_dir="$(cd "${candidate}" && pwd)"
      break
    fi
  done
fi
if [[ -z "${data_dir}" ]] || ! has_region_data "${data_dir}"; then
  echo "Could not find the 11 brain-region CSV files. Set DATA_DIR=/path/to/data." >&2
  exit 1
fi

mkdir -p "${figure_dir}" "${result_dir}"
export PROTEOMICS_FIGURE_FONT="${font_family}"

r_args=(--data-dir "${data_dir}" --output-dir "${result_dir}" --font-family "${font_family}")

"${rscript_executable}" "${root}/scripts/R/volcano_BLA_G4_vs_G2.R" "${r_args[@]}" --output-pdf "${figure_dir}/BLA_volcano_G4_vs_G2.pdf"
"${rscript_executable}" "${root}/scripts/R/go_enrichment_BLA_G4_vs_G2.R" "${r_args[@]}" --output-pdf "${figure_dir}/BLA_GO_enrichment_G4_vs_G2.pdf"
"${python_executable}" "${root}/scripts/python/brain_region_pca.py" \
  --data-dir "${data_dir}" \
  --output-dir "${result_dir}/pca" \
  --output-pdf "${figure_dir}/brain_region_pca.pdf" \
  --font-family "${font_family}"
"${rscript_executable}" "${root}/scripts/R/region_paired_heatmap_G2_G4.R" "${r_args[@]}" --output-pdf "${figure_dir}/multi_region_G2_vs_G4_heatmap.pdf"
"${rscript_executable}" "${root}/scripts/R/region_log2fc_correlation_G2_G4.R" "${r_args[@]}" --output-pdf "${figure_dir}/brain_region_log2FC_correlation.pdf" --min-pairwise-proteins 30
"${rscript_executable}" "${root}/environment/capture_session_info.R"

for name in BLA_volcano_G4_vs_G2.pdf BLA_GO_enrichment_G4_vs_G2.pdf brain_region_pca.pdf multi_region_G2_vs_G4_heatmap.pdf brain_region_log2FC_correlation.pdf; do
  test -s "${figure_dir}/${name}"
done

echo "Figure generation completed."
echo "Data directory: ${data_dir}"
echo "Figures: ${figure_dir}"
echo "Results: ${result_dir}"
