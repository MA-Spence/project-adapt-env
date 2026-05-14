#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export PYTHONPATH="${PROJECT_ROOT}/src:${PROJECT_ROOT}/external/Adapt-Env${PYTHONPATH:+:${PYTHONPATH}}"

python -u "${PROJECT_ROOT}/scripts/proteingym_phot_structural_mismatch_diagnostic.py" \
  --config "${SCRIPT_DIR}/config.yaml" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}"
