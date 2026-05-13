#!/usr/bin/env bash
#SBATCH --job-name=EXP-010RUN-030
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=/tank/data/labproj/runtime/project-adapt-env/runs/RUN-030/logs/stdout.txt
#SBATCH --error=/tank/data/labproj/runtime/project-adapt-env/runs/RUN-030/logs/stderr.txt
set -euo pipefail
export LABPROJ_RUN_ATTEMPT=1
export LABPROJ_RESUME_RUN=0
export LABPROJ_STACK_PREFIX="/tank/data/labproj/stacks/STACK-000"
export PATH="/tank/data/labproj/stacks/STACK-000/bin:${PATH}"
if ! command -v module >/dev/null 2>&1; then
  set +u
  source /etc/profile
  set -u
fi
if command -v module >/dev/null 2>&1; then
  module purge
fi
if [ -x "/tank/data/labproj/stacks/STACK-000/bin/python" ] && [ ! -x "/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin/python" ]; then
  mkdir -p "/tank/data/labproj/envs/matthew-spence/project-adapt-env"
  "/tank/data/labproj/stacks/STACK-000/bin/python" -m venv --system-site-packages "/tank/data/labproj/envs/matthew-spence/project-adapt-env"
fi
if [ -x "/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin/python" ] && [ -f "/tank/data/labproj/staging/project-adapt-env/RUN-030/project/environment/requirements.lock.txt" ]; then
  LABPROJ_OVERLAY_STAMP="/tank/data/labproj/envs/matthew-spence/project-adapt-env/.labproj-overlay.sha256"
  LABPROJ_OVERLAY_SHA=""
  if [ -f "${LABPROJ_OVERLAY_STAMP}" ]; then
    LABPROJ_OVERLAY_SHA="$(cat "${LABPROJ_OVERLAY_STAMP}")"
  fi
  if [ "${LABPROJ_OVERLAY_SHA}" != "cb34da18942b359c889eacd523eb8e584b8808a4ffe146a8e4d29c048de00859" ]; then
    "/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin/python" -m pip install --disable-pip-version-check -r "/tank/data/labproj/staging/project-adapt-env/RUN-030/project/environment/requirements.lock.txt"
    printf '%s\n' "cb34da18942b359c889eacd523eb8e584b8808a4ffe146a8e4d29c048de00859" > "${LABPROJ_OVERLAY_STAMP}"
  fi
fi
if [ -x "/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin/python" ]; then
  export LABPROJ_PROJECT_ENV_ROOT="/tank/data/labproj/envs/matthew-spence/project-adapt-env"
  export VIRTUAL_ENV="/tank/data/labproj/envs/matthew-spence/project-adapt-env"
  export PATH="/tank/data/labproj/envs/matthew-spence/project-adapt-env/bin:${PATH}"
fi
export LABPROJ_TARGET=lab-slurm
export LABPROJ_CLUSTER_NAME=homelab
export LABPROJ_EXPERIMENT_ID=EXP-010
export LABPROJ_RUN_ID=RUN-030
export LABPROJ_OUTPUT_DIR="/tank/data/labproj/runtime/project-adapt-env/runs/RUN-030/outputs"
export LABPROJ_LOG_DIR="/tank/data/labproj/runtime/project-adapt-env/runs/RUN-030/logs"
export LABPROJ_WORKDIR="${TMPDIR:-/tmp/labproj/project-adapt-env/RUN-030}"

mkdir -p "${LABPROJ_OUTPUT_DIR}" "${LABPROJ_LOG_DIR}" "${LABPROJ_WORKDIR}"

set +e
bash "/tank/data/labproj/staging/project-adapt-env/RUN-030/project/experiments/2026-05-13_EXP-010_proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/run.sh"
LABPROJ_EXIT_CODE=$?
echo "${LABPROJ_EXIT_CODE}" > "${LABPROJ_LOG_DIR}/labproj-exit-code.txt"
set -e
exit "${LABPROJ_EXIT_CODE}"
