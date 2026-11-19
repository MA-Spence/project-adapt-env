#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/local}"
DATA_ROOT="${ADAPTENV_DATA_ROOT:-/tank/data/labproj/runtime/project-adapt-env/shared/ensemble-landscape-validation-datasets}"
REPLICATE_ROOT="${ADAPTENV_REPLICATE_ROOT:-/tank/data/labproj/runtime/project-adapt-env/shared/trpb-replicate-observation-exp065}"
PRODUCTION_ROOT="${ADAPTENV_PRODUCTION_ROOT:-/tank/data/labproj/projects/project-adapt-env/data/processed/production_landscapes}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
PYTHON="${PYTHON:-/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin/python}"
set --
if [[ "${LABPROJ_QUICK:-0}" == "1" ]]; then set -- --quick; fi
"${PYTHON}" "${PROJECT_ROOT}/scripts/production_landscapes_exp066.py" \
  --config "${SCRIPT_DIR}/config.yaml" --output-dir "${OUTPUT_DIR}" \
  --data-root "${DATA_ROOT}" --replicate-data-root "${REPLICATE_ROOT}" "$@"
# Publish only a complete, hash-verified bundle; preserve any previous attempt.
"${PYTHON}" "${PROJECT_ROOT}/scripts/publish_production_landscapes.py" \
  --source "${OUTPUT_DIR}" --destination "${PRODUCTION_ROOT}/${LABPROJ_RUN_ID:-local}"
