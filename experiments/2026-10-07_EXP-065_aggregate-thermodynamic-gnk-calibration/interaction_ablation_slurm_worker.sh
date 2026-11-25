#!/usr/bin/env bash
# Frozen TrpB joint Optuna calibration, with the restored RUN-180 objective.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUNDLE="${PROJECT_ROOT}/data/processed/trpb-interaction-ablation-slurm-runtime"
PYTHON=/tank/data/labproj/runtime/project-adapt-env/envs/joint-readout-py311/bin/python
: "${LABPROJ_OUTPUT_DIR:?labproj output directory is required}"
: "${LABPROJ_LOG_DIR:?labproj log directory is required}"
if [[ "${LABPROJ_TARGET:-}" != "lab-slurm" ]]; then
  echo "This frozen worker requires the lab-slurm target." >&2
  exit 1
fi
unset PYTHONPATH PYTHONHOME ADAPTENV_DATA_ROOT
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LABPROJ_QUICK=0
mkdir -p "${LABPROJ_LOG_DIR}"
"${PYTHON}" "${SCRIPT_DIR}/submit_interaction_ablation_slurm.py" \
  --worker-preflight "${BUNDLE}" --log-dir "${LABPROJ_LOG_DIR}"
echo "Starting matched TrpB r0/full shared-stop calibration; RUN=${LABPROJ_RUN_ID:-unknown}; job=${SLURM_JOB_ID:-unknown}"
set +e
"${PYTHON}" "${BUNDLE}/source/scripts/interaction_ablation_exp065.py" \
  --config "${BUNDLE}/source/experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/interaction_ablation.yaml" \
  --output-dir "${LABPROJ_OUTPUT_DIR}" \
  --baseline-root "${BUNDLE}/baseline/RUN-186"
STATUS=$?
set -e
printf '%s\n' "${STATUS}" > "${LABPROJ_LOG_DIR}/interaction-ablation-scientific-exit-code.txt"
exit "${STATUS}"
