#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/demo}"
WORK_DIR="${LABPROJ_WORKDIR:-${SCRIPT_DIR}/work/demo}"
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"
python "${PROJECT_ROOT}/scripts/proteingym_latent_observation_calibration.py" \
  --config "${SCRIPT_DIR}/config.yaml" \
  --output-dir "${OUTPUT_DIR}" \
  --project-root "${PROJECT_ROOT}"
