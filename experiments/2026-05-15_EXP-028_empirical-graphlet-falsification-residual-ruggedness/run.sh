#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
CONFIG_PATH="${SCRIPT_DIR}/config.yaml"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

export PYTHONUNBUFFERED=1
export PYTHONPATH="${PROJECT_ROOT}/src:${PROJECT_ROOT}/external/Adapt-Env:${PROJECT_ROOT}/scripts${PYTHONPATH:+:${PYTHONPATH}}"

QUICK_ARGS=()
if [[ "${LABPROJ_QUICK:-0}" == "1" ]]; then
  QUICK_ARGS+=(--quick)
fi

"${PYTHON:-python3}" -u "${PROJECT_ROOT}/scripts/empirical_graphlet_falsification_suite.py" \
  --config "${CONFIG_PATH}" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}" \
  "${QUICK_ARGS[@]}"
