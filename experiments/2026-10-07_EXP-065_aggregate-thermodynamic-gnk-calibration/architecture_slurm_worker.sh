#!/usr/bin/env bash
# Transport-only worker. Scientific source, settings and inputs remain RUN-182's.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUNDLE="${PROJECT_ROOT}/data/processed/trpb-architecture-slurm-runtime"
PYTHON=/tank/data/labproj/runtime/project-adapt-env/envs/architecture-0.7.0rc1-py311/bin/python
: "${LABPROJ_OUTPUT_DIR:?labproj output directory is required}"
: "${LABPROJ_LOG_DIR:?labproj log directory is required}"
if [[ "${LABPROJ_TARGET:-}" != "lab-slurm" ]]; then
  echo "This frozen worker requires the lab-slurm target." >&2
  exit 1
fi
unset PYTHONPATH PYTHONHOME
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LABPROJ_QUICK=0
mkdir -p "${LABPROJ_LOG_DIR}"
"${PYTHON}" "${SCRIPT_DIR}/submit_architecture_slurm.py" \
  --worker-preflight "${BUNDLE}" --log-dir "${LABPROJ_LOG_DIR}"
echo "Starting unchanged RUN-182 calibration on Linux; RUN=${LABPROJ_RUN_ID:-unknown}; job=${SLURM_JOB_ID:-unknown}"
set +e
"${PYTHON}" "${BUNDLE}/source/scripts/calibrate_architecture.py" \
  --input "${BUNDLE}/inputs/trpb-architecture-calibration-inputs/calibration.npz" \
  --confirmation "${BUNDLE}/inputs/trpb-architecture-calibration-inputs/confirmation.npz" \
  --config "${BUNDLE}/source/experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/architecture_calibration.json" \
  --output "${LABPROJ_OUTPUT_DIR}"
STATUS=$?
set -e
printf '%s\n' "${STATUS}" > "${LABPROJ_LOG_DIR}/architecture-scientific-exit-code.txt"
exit "${STATUS}"
