#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

export PYTHONPATH="${PROJECT_ROOT}/external/Adapt-Env:${PROJECT_ROOT}/scripts:${PYTHONPATH:-}"

python "${PROJECT_ROOT}/scripts/tem1_dca_clique_haplotype_epistasis.py" \
  --config "${SCRIPT_DIR}/config.yaml" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}"
