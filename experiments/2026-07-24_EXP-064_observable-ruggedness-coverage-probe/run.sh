#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
CONFIG_PATH="${SCRIPT_DIR}/config.yaml"
RUNNER="${PROJECT_ROOT}/scripts/observable_ruggedness_coverage_probe_exp064.py"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export LABPROJ_WORKERS="${LABPROJ_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
export PYTHONPATH="${PROJECT_ROOT}/src:${PROJECT_ROOT}/external/Adapt-Env:${PROJECT_ROOT}/scripts${PYTHONPATH:+:${PYTHONPATH}}"

DATA_DIR="${PROJECT_ROOT}/data/processed/ensemble-landscape-validation-datasets"
DATA_DVC="${PROJECT_ROOT}/data/processed/ensemble-landscape-validation-datasets.dvc"
SHARED_DATA_DIR="/tank/data/labproj/runtime/project-adapt-env/shared/ensemble-landscape-validation-datasets"
if [[ ! -d "${DATA_DIR}/datasets" && -d "${SHARED_DATA_DIR}/datasets" ]]; then
  export STAGE4_PROCESSED_ROOT="${SHARED_DATA_DIR}"
elif [[ ! -d "${DATA_DIR}/datasets" && -f "${DATA_DVC}" ]]; then
  (cd "${PROJECT_ROOT}" && dvc pull data/processed/ensemble-landscape-validation-datasets.dvc)
fi

PYTHON_BIN="${PYTHON:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${PROJECT_ROOT}/.venv/bin/python"
  else
    PYTHON_BIN="python3"
  fi
fi

QUICK_ARGS=()
if [[ "${LABPROJ_QUICK:-0}" == "1" ]]; then
  QUICK_ARGS+=(--quick)
fi

"${PYTHON_BIN}" -u "${RUNNER}" \
  --config "${CONFIG_PATH}" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}" \
  "${QUICK_ARGS[@]}"
