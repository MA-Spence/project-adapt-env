#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
CONFIG_PATH="${SCRIPT_DIR}/config.yaml"
SOURCE_DVC="${PROJECT_ROOT}/data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079.dvc"
SOURCE_DIR="${PROJECT_ROOT}/data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export PYTHONPATH="${PROJECT_ROOT}/src:${PROJECT_ROOT}/external/Adapt-Env:${PROJECT_ROOT}/scripts${PYTHONPATH:+:${PYTHONPATH}}"

if [[ ! -d "${SOURCE_DIR}" && -f "${SOURCE_DVC}" ]] && command -v dvc >/dev/null 2>&1; then
  dvc pull "${SOURCE_DVC}"
fi

QUICK_ARGS=()
if [[ "${LABPROJ_QUICK:-0}" == "1" ]]; then
  QUICK_ARGS+=(--quick)
fi

"${PYTHON:-python3}" -u "${PROJECT_ROOT}/scripts/phot_replicated_split_graphlet_diagnostic.py" \
  --config "${CONFIG_PATH}" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}" \
  "${QUICK_ARGS[@]}"
