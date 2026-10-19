#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
args=("${root}/run_all.py" --data-dir "${DATA_DIR:?Set DATA_DIR to the directory of 12 original CSV reports}" --rscript "${RSCRIPT:-Rscript}" --figures "${FIGURES:-all}" --environment "${ENVIRONMENT:-manuscript}")
if [[ -n "${METADATA:-}" ]]; then args+=(--metadata "${METADATA}"); fi
if [[ -n "${R_LIBS_USER:-}" ]]; then args+=(--r-library "${R_LIBS_USER}"); fi
if [[ -n "${OUTPUT_ROOT:-}" ]]; then args+=(--output-root "${OUTPUT_ROOT}"); fi
"${PYTHON:-python3}" "${args[@]}" "$@"
